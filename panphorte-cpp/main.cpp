#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <chrono> // For performance timing

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <fcntl.h> // For redirecting stdout/stderr

#include "PangenomeGraph.hpp"
#include "BubbleProcessor.hpp"

namespace fs = std::filesystem;

// -----------------------------------------------------------------------------
// CLI & System Helpers
// -----------------------------------------------------------------------------

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

    if (pid == -1) {
        // Fork failed
        std::cerr << "[fatal] Failed to fork process for BubbleGun." << std::endl;
        return -1;
    } else if (pid == 0) {
        // --- CHILD PROCESS ---
        
        // 2. Redirect stdout/stderr (Optional, matches your "> /dev/null")
        // Opening /dev/null
        int dev_null = open("/dev/null", O_WRONLY);
        if (dev_null != -1) {
            dup2(dev_null, STDOUT_FILENO); // Redirect stdout
            dup2(dev_null, STDERR_FILENO); // Redirect stderr
            close(dev_null);
        }

        // 3. Prepare Arguments safely
        // The array must be NULL terminated.
        // const_cast is safe here because execvp doesn't modify args 
        // (though standard signature is char* const[])
        std::vector<const char*> args;
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
        execvp("BubbleGun", const_cast<char* const*>(args.data()));

        // execvp failed (e.g., tool not found)
        std::cerr << "[fatal] Could not find or execute 'BubbleGun'. Is it in your PATH?" << std::endl;
        _exit(127); // to avoid flushing parent's buffers
    } else {
        // --- PARENT PROCESS ---
        
        // 5. Wait for child to finish
        int status;
        waitpid(pid, &status, 0);

        if (WIFEXITED(status)) {
            int exit_code = WEXITSTATUS(status);
            if (exit_code == 0) {
                std::cout << "[exec] BubbleGun finished successfully." << std::endl;
                return 0;
            } else {
                std::cerr << "[error] BubbleGun failed with exit code: " << exit_code << std::endl;
                return exit_code;
            }
        } else {
            std::cerr << "[error] BubbleGun terminated abnormally." << std::endl;
            return -1;
        }
    }
}

int main(int argc, char **argv)
{
    // 1. Parse Arguments
    std::string gfa_path;
    std::string output_dir;
    if (!parse_args(argc, argv, gfa_path, output_dir)) {
        usage(argv[0]);
        return 1;
    }

    // 2. Setup Paths
    try {
        if (!fs::exists(output_dir)) {
            fs::create_directories(output_dir);
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Error creating output directory: " << e.what() << std::endl;
        return 1;
    }

    fs::path in_p(gfa_path);
    std::string stem = in_p.stem().string();
    fs::path out_gfa_path = fs::path(output_dir) / (stem + "_mod.gfa");
    
    // BubbleGun intermediate files
    std::string bubble_json = (fs::path(output_dir) / "bubbles.json").string();
    std::string bubble_fasta = (fs::path(output_dir) / "bubbles.fasta").string();

    // 3. Initialize and Load Graph
    PangenomeGraph graph;

    std::cout << "[info] Loading GFA: " << gfa_path << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();

    try {
        graph.load_from_gfa(gfa_path);
    } 
    catch (const std::exception& e) {
        std::cerr << "[fatal] Failed to parse GFA: " << e.what() << std::endl;
        return 1;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    // 4. Print Statistics (Validation)
    std::cout << "[info] Loaded successfully in " << elapsed.count() << "s\n";
    // std::cout << "       - Header lines: " << graph.header.size() << "\n";
    std::cout << "       - Nodes: " << graph.get_node_count() << "\n";
    std::cout << "       - Links: " << graph.links.size() << "\n";
    std::cout << "       - Walks: " << graph.walks.size() << "\n";

    // 5. Run BubbleGun (External Tool)
    if (run_bubblegun(gfa_path, bubble_json, bubble_fasta) != 0) {
        return 1;
    }

    // 6. Parse and Save Information
    std::vector<BubbleContext> bubble_tasks;
    try {
        bubble_tasks = parse_bubblegun_json(bubble_json, graph);
    } catch (const std::exception& e) {
        std::cerr << "[fatal] JSON Parsing failed: " << e.what() << std::endl;
        return 1;
    }

    

    return 0;
}