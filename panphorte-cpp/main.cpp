#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <chrono> // For performance timing
#include <omp.h>

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <fcntl.h> // For redirecting stdout/stderr

#include "PangenomeGraph.hpp"
#include "BubbleProcessor.hpp"

namespace fs = std::filesystem;

static void usage(const char *prog)
{
    std::cerr << "Usage: " << prog << " -i <graph.gfa> [-o <output_dir>]\n"
              << "  -i, --input        Input GFA file (required)\n"
              << "  -o, --output_dir   Output directory for the modified GFA (default: .)\n";
}

static bool parse_args(int argc, char **argv, std::string &gfa_path, std::string &output_dir)
{
    output_dir = ".";
    for (int i = 1; i < argc; ++i)
    {
        std::string a = argv[i];
        if ((a == "-i" || a == "--input") && i + 1 < argc)
        {
            gfa_path = argv[++i];
            if (gfa_path.size() < 4 || gfa_path.substr(gfa_path.size() - 4) != ".gfa")
            {
                std::cerr << "Error: input file must have .gfa extension\n";
                return false;
            }
        }
        else if ((a == "-o" || a == "--output_dir") && i + 1 < argc)
        {
            output_dir = argv[++i];
        }
    }
    return !gfa_path.empty();
}

int run_bubblegun(const std::string &gfa_file,
                  const std::string &out_json,
                  const std::string &fasta_out)
{
    std::cout << "[exec] Running BubbleGun..." << std::endl;

    // 1. Fork the process
    pid_t pid = fork();

    if (pid == -1)
    {
        // Fork failed
        std::cerr << "[fatal] Failed to fork process for BubbleGun." << std::endl;
        return -1;
    }
    else if (pid == 0)
    {
        // --- CHILD PROCESS ---

        // 2. Redirect stdout/stderr (Optional, matches your "> /dev/null")
        // Opening /dev/null
        int dev_null = open("/dev/null", O_WRONLY);
        if (dev_null != -1)
        {
            dup2(dev_null, STDOUT_FILENO); // Redirect stdout
            dup2(dev_null, STDERR_FILENO); // Redirect stderr
            close(dev_null);
        }

        // 3. Prepare Arguments safely
        // The array must be NULL terminated.
        // const_cast is safe here because execvp doesn't modify args
        // (though standard signature is char* const[])
        std::vector<const char *> args;
        args.push_back("BubbleGun");
        args.push_back("-g");
        args.push_back(gfa_file.c_str());
        args.push_back("bchains");
        args.push_back("--bubble_json");
        args.push_back(out_json.c_str());
        args.push_back("--fasta");
        args.push_back(fasta_out.c_str());
        args.push_back(nullptr); // Terminator

        // 4. Execute
        // execvp searches PATH for "BubbleGun"
        execvp("BubbleGun", const_cast<char *const *>(args.data()));

        // execvp failed (e.g., tool not found)
        std::cerr << "[fatal] Could not find or execute 'BubbleGun'. Is it in your PATH?" << std::endl;
        _exit(127); // to avoid flushing parent's buffers
    }
    else
    {
        // --- PARENT PROCESS ---

        // 5. Wait for child to finish
        int status;
        waitpid(pid, &status, 0);

        if (WIFEXITED(status))
        {
            int exit_code = WEXITSTATUS(status);
            if (exit_code == 0)
            {
                std::cout << "[exec] BubbleGun finished successfully." << std::endl;
                return 0;
            }
            else
            {
                std::cerr << "[error] BubbleGun failed with exit code: " << exit_code << std::endl;
                return exit_code;
            }
        }
        else
        {
            std::cerr << "[error] BubbleGun terminated abnormally." << std::endl;
            return -1;
        }
    }
}

int main(int argc, char **argv)
{
    std::string gfa_path;
    std::string output_dir;
    if (!parse_args(argc, argv, gfa_path, output_dir))
    {
        usage(argv[0]);
        return 1;
    }

    try
    {
        if (!fs::exists(output_dir))
            fs::create_directories(output_dir);
    }
    catch (const fs::filesystem_error &e)
    {
        std::cerr << "[fatal] Error creating output directory: " << e.what() << std::endl;
        return 1;
    }

    fs::path in_p(gfa_path);
    std::string stem = in_p.stem().string();
    fs::path out_gfa_path = fs::path(output_dir) / (stem + "_mod_new.gfa");
    std::string bubble_json = (fs::path(output_dir) / "bubbles.json").string();
    std::string bubble_fasta = (fs::path(output_dir) / "bubbles.fasta").string();

    PangenomeGraph graph;
    std::cout << "[info] Loading GFA: " << gfa_path << std::endl;
    auto start_load = std::chrono::high_resolution_clock::now();

    try
    {
        graph.load_from_gfa(gfa_path);
    }
    catch (const std::exception &e)
    {
        std::cerr << "[fatal] Failed to parse GFA: " << e.what() << std::endl;
        return 1;
    }

    // Allows O(1) haplotype extraction later
    graph.build_reverse_index();

    auto end_load = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_load = end_load - start_load;

    // Statistics
    std::cout << "[info] Loaded in " << elapsed_load.count() << "s\n";
    std::cout << "       - Nodes: " << graph.get_node_count() << "\n";
    std::cout << "       - Walks: " << graph.walks.size() << "\n";
    std::cout << "       - Links: " << graph.links.size() << "\n";

    // Run BubbleGun
    if (run_bubblegun(gfa_path, bubble_json, bubble_fasta) != 0)
    {
        return 1;
    }

    // Parse BubbleGun JSON
    std::vector<BubbleContext> bubble_tasks;
    try
    {
        bubble_tasks = parse_bubblegun_json(bubble_json, graph);
    }
    catch (const std::exception &e)
    {
        std::cerr << "[fatal] JSON Parsing failed: " << e.what() << std::endl;
        return 1;
    }

    // Parallel Analysis (Step 7)
    // We separate analysis (Read-Only) from Modification (Write)
    std::vector<ModificationTask> modifications;

    std::cout << "[info] Analyzing " << bubble_tasks.size() << " bubbles for repetitions...\n";
    auto start_analysis = std::chrono::high_resolution_clock::now();

    int found_count = 0;

#pragma omp parallel for schedule(dynamic) reduction(+ : found_count)
    for (size_t i = 0; i < bubble_tasks.size(); ++i)
    {
        const auto &bubble = bubble_tasks[i];

        // std::cout << "[debug] Processing Bubble " << i
        //           << " (Start: " << graph.get_node_name(bubble.start_node)
        //           << ", End: " << graph.get_node_name(bubble.end_node) << ")\n";

        std::vector<Haplotype> haps = extract_bubble_haplotypes(bubble, graph);

        // print haps for debugging
//         for (const auto &h : haps)
//         {
// #pragma omp critical
//             {
//                 std::cout << "[debug] Bubble " << i << " Haplotype: " << h.sequence
//                           << " (Nodes: ";
//                 for (auto nid : h.path_nodes)
//                 {
//                     std::cout << graph.get_node_name(nid) << " ";
//                 }
//                 std::cout << ")\n";
//             }
//         }

        auto task_opt = find_common_repetition(haps, bubble.start_node, bubble.end_node);

        if (task_opt.has_value())
        {
// Critical section to safely push to shared vector
#pragma omp critical
            {
                modifications.push_back(task_opt.value());
            }
            found_count++;
        }
    }

    auto end_analysis = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_analysis = end_analysis - start_analysis;
    std::cout << "[info] Analysis done in " << elapsed_analysis.count() << "s\n";
    std::cout << "[info] Found " << found_count << " optimization targets.\n";

    // Apply Modifications
    std::cout << "[info] Applying topology optimizations to " << modifications.size() << " nodes...\n";

    int applied_count = 0;
    for (const auto &mod : modifications)
    {
        // Safety check: Ensure the node hasn't been deleted by a previous overlapping modification
        // (Rare)
        if (!graph.is_node_deleted(mod.bubble_start) && !graph.is_node_deleted(mod.bubble_end))
        {
            graph.apply_bubble_optimization(mod);
            applied_count++;
        }
    }
    std::cout << "[info] Successfully optimized " << applied_count << " regions.\n";

    // Save Result
    std::cout << "[info] Saving GFA to " << out_gfa_path.string() << "...\n";
    graph.save_to_gfa(out_gfa_path.string());
    std::cout << "[info] Done.\n";

    return 0;
}