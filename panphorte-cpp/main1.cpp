#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <cstdio> // popen, pclose, fgets
#include <array>
#include <cerrno>
#include <sys/wait.h> // WIFEXITED, WEXITSTATUS (POSIX)
#include <regex>

#include <sstream>
#include <cstdlib>
#include <stdexcept>
#include "json.hpp"

#include <optional>
#include <algorithm>
#include <random>

#include <chrono>

using json = nlohmann::json;
using namespace std;

struct Link
{
    string source;
    string source_orient;
    string target;
    string target_orient;
    string overlap;
};

struct Walk
{
    std::string sample;
    int hap; // forse meglio string
    std::string seqid;
    std::string start;
    std::string end;
    std::vector<std::pair<std::string, char>> segments; // (segId, orient) con orient in {'+','-'}
};

struct GFAData
{
    vector<string> header;
    unordered_map<string, string> nodes;         // id -> sequence
    unordered_map<string, vector<string>> paths; // path_id -> segments (senza orientamento)
    vector<Link> links;

    unordered_map<string, Walk> walks; // walk_id -> Walk
};

// ---- Utility functions for Paths ----
static inline void split_tab(const string &s, vector<string> &out)
{
    out.clear();
    stringstream ss(s);
    string item;
    while (getline(ss, item, '\t'))
        out.push_back(item);
}

static inline void split_comma(const string &s, vector<string> &out)
{
    out.clear();
    stringstream ss(s);
    string item;
    while (getline(ss, item, ','))
        out.push_back(item);
}

static inline bool ends_with_plus_or_minus(const string &s)
{
    if (s.empty())
        return false;
    char c = s.back();
    return c == '+' || c == '-';
}

// ---- Utility functions for Walks ----
static inline string trim(const string &s)
{
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == string::npos)
        return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

// Parse W field tipo ">s11<s12>s13" in lista di (segId, orient) con orient in {'+','-'}
static inline vector<pair<string, char>> parse_walk_string(const string &walk_str)
{
    static const regex token_re(R"(([<>])([^<>\s,]+))");
    vector<pair<string, char>> result;
    smatch m;
    string s = trim(walk_str);
    auto it = s.cbegin();
    auto ed = s.cend();

    while (regex_search(it, ed, m, token_re))
    {
        const string arrow = m[1].str();
        const string seg = m[2].str();
        const char orient = (arrow == ">") ? '+' : '-';
        result.emplace_back(seg, orient);
        it = m.suffix().first;
    }
    return result;
}

static std::string format_walk_segments(const std::vector<std::pair<std::string, char>> &segs)
{
    std::string out;
    out.reserve(segs.size() * 8);
    for (const auto &[seg, orient] : segs)
    {
        out.push_back(orient == '+' ? '>' : '<');
        out += seg;
    }
    return out;
}

// Helpers per il parsing del json
// Restituisce true se segments contiene un tuple (seg_id, _)
// Not used anymore!!!
// static bool contains_seg(const std::vector<std::pair<std::string, char>> &segments,
//                          const std::string &seg_id)
// {
//     for (const auto &[sid, orient] : segments)
//     {
//         if (sid == seg_id)
//             return true;
//     }
//     return false;
// }

// Restituisce l’indice del primo segmento che ha seg_id == sid
// Se non lo trova, solleva un’eccezione (come il tuo ValueError)
// int find_index_by_seg(const std::vector<std::pair<std::string,char>>& segments,
//                       const std::string& seg_id)
// {
//     for (size_t i = 0; i < segments.size(); ++i) {
//         if (segments[i].first == seg_id)
//             return static_cast<int>(i);
//     }
//     throw std::runtime_error(seg_id + " not found in walk");
// }

static int find_index_by_seg(const std::vector<std::pair<std::string, char>> &segments,
                             const std::string &seg_id)
{
    for (size_t i = 0; i < segments.size(); ++i)
        if (segments[i].first == seg_id)
            return static_cast<int>(i);
    return -1; // not found
}

// with fixed size random ids (at most, 10 characters after prefix)
static std::string new_random_id(const std::unordered_map<std::string, std::string> &existing, const std::string &prefix)
{
    static std::mt19937_64 rng{std::random_device{}()};
    for (;;)
    {
        uint64_t x_full = rng();

        // 1. Reduce entropy to fit in 32 bits (4,294,967,295 max value).
        // This guarantees the resulting decimal string will be at most 10 digits.
        uint32_t x_short = static_cast<uint32_t>(x_full);

        // 2. Convert the shortened number to a string.
        std::string x_str = std::to_string(x_short);

        std::string id = prefix + "_" + x_str;

        if (!existing.count(id))
            return id;
    }
}

// genera un ID unico col prefisso dato, evitando collisioni in `nodes` ('existing')
// static std::string new_random_id(const std::unordered_map<std::string, std::string> &existing, const std::string &prefix)
// {
//     // TODO: accorciare lunghezza dell'ID
//     static std::mt19937_64 rng{std::random_device{}()};
//     for (;;)
//     {
//         uint64_t x = rng();
//         std::string id = prefix + "_" + std::to_string(x);
//         if (!existing.count(id))
//             return id;
//     }
// }

// test function of new_random_id for debugging with incremental ids
// int test_id = 1;
// static std::string new_random_id(const std::unordered_map<std::string, std::string> &existing, const std::string &prefix)
// {
//     for (;;)
//     {
//         std::string id = prefix + "_" + std::to_string(test_id++);
//         if (!existing.count(id))
//             return id;
//     }
// }

// Ritorna: vettore di tuple (motif, times, start_pos)
using Repetition = std::tuple<std::string, int, int>;
static std::vector<Repetition> regex_find_repetitions(const std::string &seq)
{
    std::vector<Repetition> reps;
    if (seq.empty())
        return reps;

    const int repetition_length = 1; // mirror Python's repetition_length (adjust if needed)
    try
    {
        std::regex pat(R"((.+?)\1+)");
        auto it = std::sregex_iterator(seq.begin(), seq.end(), pat);
        auto end = std::sregex_iterator();
        for (; it != end; ++it)
        {
            const std::smatch &m = *it;
            std::string motif = m[1].str();
            if ((int)motif.size() > repetition_length)
            {
                int total_len = static_cast<int>(m.str().size()); // full matched substring length
                int times = total_len / static_cast<int>(motif.size());
                int position = static_cast<int>(m.position(0)); // start position of full match
                reps.emplace_back(motif, times, position);
            }
        }
    }
    catch (const std::regex_error &)
    {
        // On regex errors, return empty vector
    }
    return reps;
}

static void usage(const char *prog)
{
    cerr << "Usage: " << prog << " -i <graph.gfa> [-o <output_dir>]\n"
         << "  -i, --input       Input GFA file (required)\n"
         << "  -o, --output_dir  Output directory for the modified GFA (default: .)\n";
}

static bool parse_args(int argc, char **argv, string &gfa_path, string &output_dir)
{
    output_dir = "."; // default
    for (int i = 1; i < argc; ++i)
    {
        string a = argv[i];
        if ((a == "-i" || a == "--input") && i + 1 < argc)
        {
            gfa_path = argv[++i];
            if (gfa_path.size() < 4 || gfa_path.substr(gfa_path.size() - 4) != ".gfa")
            {
                cerr << "Error: input file must have .gfa extension\n";
                return false;
            }
        }
        else if ((a == "-o" || a == "--output_dir") && i + 1 < argc)
        {
            output_dir = argv[++i];
        }
        else
        {
            cerr << "Unknown or malformed option: " << a << "\n";
            return false;
        }
    }
    return !gfa_path.empty();
}

int run_bubblegun(const std::string &bubblegun_bin,
                  const std::string &gfa_file,
                  const std::string &out_json,
                  const std::string &fasta,
                  const std::string &extra = "",
                  const std::string &log_stdout = "/tmp/bubblegun_stdout.log",
                  const std::string &log_stderr = "/tmp/bubblegun_stderr.log")
{
    // Python: ["BubbleGun", "-g", gfa_file, "bchains", "--bubble_json", out_json, "--fasta", fasta]
    std::string cmd = bubblegun_bin +
                      " -g " + gfa_file +
                      " bchains " +
                      " --bubble_json " + out_json +
                      " --fasta " + fasta +
                      (extra.empty() ? "" : (" " + extra)) +
                      " 1> " + log_stdout + " 2> " + log_stderr;

    int rc = std::system(cmd.c_str());
    if (rc == 0)
    {
        std::cout << "Command executed with success!\n";
        std::cout << "See stdout: " << log_stdout << "\n";
    }
    else
    {
        std::cerr << "Error with the execution of the command. rc=" << rc << "\n";
        std::cerr << "See stderr: " << log_stderr << "\n";
    }
    return rc;
}

GFAData read_gfa(const string &file_gfa)
{
    GFAData out;
    std::ifstream in(file_gfa);
    if (!in)
    {
        throw runtime_error("Impossibile aprire il file GFA: " + file_gfa);
    }

    std::string line;
    std::vector<std::string> fields;
    size_t walk_auto_id = 0; // per generare walk_id come in Python
    while (getline(in, line))
    {
        if (line.empty())
            continue;

        // Ignora eventuali righe di commento complete
        if (line[0] == '#')
            continue;

        char tipo_record = line[0];
        switch (tipo_record)
        {
        case 'H':
        {
            // Mantieni la riga intera
            out.header.push_back(line);
            break;
        }
        case 'S':
        {
            // Formato minimo atteso: S <sid> <sequence> ...
            split_tab(line, fields);
            if (fields.size() >= 3)
            {
                const std::string &id_node = fields[1];
                const std::string &sequence = fields[2];
                out.nodes[id_node] = sequence;
            }
            break;
        }
        case 'P':
        {
            // Formato minimo atteso: P <pid> <segment(+/-),...> <overlaps or *> ...
            split_tab(line, fields);
            if (fields.size() >= 3)
            {
                const std::string &path_id = fields[1];
                std::vector<std::string> segs_raw;
                split_comma(fields[2], segs_raw);

                std::vector<std::string> segs_clean;
                segs_clean.reserve(segs_raw.size());
                for (const auto &s : segs_raw)
                {
                    if (!s.empty() && ends_with_plus_or_minus(s))
                    {
                        segs_clean.emplace_back(s.substr(0, s.size() - 1));
                    }
                    else
                    {
                        segs_clean.emplace_back(s);
                    }
                }
                out.paths[path_id] = std::move(segs_clean);
            }
            break;
        }
        case 'L':
        {
            // Formato minimo atteso: L <from> <from_orient> <to> <to_orient> <overlap> ...
            split_tab(line, fields);
            if (fields.size() >= 6)
            { // nel codice python è 5, ma servono almeno 6 campi
                Link lk;
                lk.source = fields[1];
                lk.source_orient = fields[2];
                lk.target = fields[3];
                lk.target_orient = fields[4];
                lk.overlap = fields[5];
                out.links.push_back(std::move(lk));
            }
            break;
        }
        case 'W':
        {
            // W  SampleId  HapIndex  SeqId  SeqStart  SeqEnd  Walk
            split_tab(line, fields);
            if (fields.size() >= 7)
            {
                const std::string &sample = fields[1];
                const std::string &hap_str = fields[2];
                const std::string &seqid = fields[3];
                const std::string &seqstart = fields[4];
                const std::string &seqend = fields[5];
                const std::string &walk_field = fields[6];

                auto segments = parse_walk_string(walk_field);

                // costruzione del walk_id, composto da:
                // walk_id = f"W_{walk_auto_id}:{sample}|{hap}|{seqid}|{seqstart}|{seqend}"
                // NOTA: questo id potrebbe essere troppo "verboso" e di conseguenza è uno spreco di memoria
                // (se le walk sono tante)
                std::string walk_id;
                walk_id.reserve(32 + sample.size() + seqid.size() + seqstart.size() + seqend.size());
                walk_id += "W_";
                walk_id += to_string(walk_auto_id);
                walk_id += ":";
                walk_id += sample;
                walk_id += "|";
                walk_id += hap_str;
                walk_id += "|";
                walk_id += seqid;
                walk_id += "|";
                walk_id += seqstart;
                walk_id += "|";
                walk_id += seqend;

                ++walk_auto_id;

                Walk w;
                w.sample = sample;
                w.hap = stoi(hap_str);
                w.seqid = seqid;
                w.start = seqstart;
                w.end = seqend;
                segments.shrink_to_fit();
                w.segments = std::move(segments);

                out.walks.emplace(std::move(walk_id), std::move(w));
            }
            break;
        }
        default:
            // Altri record GFA ignorati
            break;
        }
    }

    return out;
}

// Funzione helper per decidere se invertire l'iterazione
// Ritorna TRUE se dobbiamo iterare al contrario (dall'ultima alla prima)
bool should_iterate_backwards(const json& chain, const GFAData& gfa_data) {
    if (chain["bubbles"].size() < 2) return false;

    // 1. Troviamo una walk di riferimento (qualsiasi walk che passi per la catena)
    //    Ci serve solo per capire dov'è l'inizio biologico.
    const std::string* ref_walk_ptr = nullptr;
    std::string start_node = chain["ends"][0];
    
    // Cerca velocemente una walk che contiene il nodo start della catena
    if (gfa_data.nodes.empty()) return false;
    
    // Usa l'indice inverso se disponibile, altrimenti scan rapido
    // Qui assumiamo scan rapido per semplicità, ci si ferma alla prima trovata
    for (const auto& [w_id, walk] : gfa_data.walks) {
        for (const auto& seg : walk.segments) {
            if (seg.first == start_node) {
                ref_walk_ptr = &w_id;
                break;
            }
        }
        if (ref_walk_ptr) break;
    }
    
    if (!ref_walk_ptr) return false; // Non possiamo determinare, default forward

    // 2. Analizziamo la PRIMA bolla del JSON
    const json& first_bubble = chain["bubbles"][0];
    std::string b_end1 = first_bubble["ends"][0];
    std::string b_end2 = first_bubble["ends"][1];

    // 3. Analizziamo l'ULTIMA bolla del JSON
    const json& last_bubble = chain["bubbles"].back();
    std::string l_end1 = last_bubble["ends"][0];
    std::string l_end2 = last_bubble["ends"][1];

    // 4. Verifichiamo la posizione genomica sulla walk di riferimento
    const auto& segs = gfa_data.walks.at(*ref_walk_ptr).segments;
    
    long long pos_first_bubble = 999999999;
    long long pos_last_bubble = 999999999;

    for (size_t i = 0; i < segs.size(); ++i) {
        if (pos_first_bubble == 999999999 && (segs[i].first == b_end1 || segs[i].first == b_end2)) {
            pos_first_bubble = (long long)i;
        }
        if (pos_last_bubble == 999999999 && (segs[i].first == l_end1 || segs[i].first == l_end2)) {
            pos_last_bubble = (long long)i;
        }
        if (pos_first_bubble != 999999999 && pos_last_bubble != 999999999) break;
    }

    // SE la prima bolla nel JSON appare DOPO l'ultima bolla nel genoma
    // ALLORA il JSON è "al contrario" -> Dobbiamo iterare Backwards.
    return pos_first_bubble > pos_last_bubble;
}

void process_bubblegun_output(const string &file_path, GFAData &gfa_data)
{
    // --- TIMING VARIABLES ACCUMULATORS ---
    long long total_time_haplo_scan = 0;
    long long total_time_rep_find = 0;
    long long total_time_graph_mod = 0;
    // -------------------------------------

    // --- TIMING SECTION: JSON LOAD ---
    auto start_json = std::chrono::high_resolution_clock::now();

    std::ifstream in(file_path);
    if (!in)
        throw runtime_error("Cannot open JSON: " + file_path);
    json data;
    in >> data;

    auto end_json = std::chrono::high_resolution_clock::now();
    auto duration_json = std::chrono::duration_cast<std::chrono::microseconds>(end_json - start_json);
    std::cerr << "[info] JSON Parsing executed in " << duration_json.count() / 1e6 << " seconds.\n";
    // ---------------------------------

    // Mappa dei cursori: per ogni walk, ricorda l'ultimo indice visitato (inizialmente 0)
    std::unordered_map<std::string, size_t> walk_cursors;
    for (const auto& [w_id, _] : gfa_data.walks) {
        walk_cursors[w_id] = 0;
    }

    int innn = 0;
    for (auto it = data.begin(); it != data.end(); ++it)
    {
        const json &chain = it.value();
        if (!chain.contains("bubbles") || !chain["bubbles"].is_array())
            continue;

        bool backwards = true;

        // bool backwards = should_iterate_backwards(chain, gfa_data);

        size_t num_bubbles = chain["bubbles"].size();
        
        for (size_t i = 0; i < num_bubbles; ++i)
        {
            // Se backwards, prendiamo l'indice partendo dal fondo
            size_t idx = backwards ? (num_bubbles - 1 - i) : i;
            
            const json &bubble = chain["bubbles"][idx];

            ++innn;
            // cout << "Processing bubble " << bubble["id"] << "...\n";
            if (!bubble.contains("ends") || !bubble["ends"].is_array() || bubble["ends"].size() < 2)
                continue;

            std::string start_node = bubble["ends"][0].get<std::string>();
            std::string final_node = bubble["ends"][1].get<std::string>();

            // Identify correct order of ends by scanning walks
            bool ok = false;
            for (auto &kv : gfa_data.walks)
            {
                const auto &segs = kv.second.segments;
                for (auto &pr : segs)
                {
                    const auto &sid = pr.first;
                    if (sid == start_node)
                    {
                        ok = true;
                        break;
                    }
                    if (sid == final_node)
                    {
                        std::swap(final_node, start_node);
                        ok = true;
                        break;
                    }
                }
                if (ok)
                    break;
            }

            if (!bubble.contains("inside") || !bubble["inside"].is_array())
                continue;

            std::vector<std::vector<std::string>> haplotypes;

            // --- TIMING SECTION START: Haplotype Scanning ---
            auto start_scan = std::chrono::high_resolution_clock::now();

            for (const auto& kv : gfa_data.walks) 
            {
                const std::string& w_id = kv.first;
                // Otteniamo riferimento alla walk e al suo cursore
                const auto& walk = kv.second;
                size_t& current_cursor = walk_cursors[w_id]; // Riferimento per aggiornarlo dopo

                // print current cursor and w_id for debugging
                // std::cout << "Current cursor for walk " << w_id << ": " << current_cursor << "\n";

                // Se il cursore è già alla fine della walk, saltiamo
                if (current_cursor >= walk.segments.size()) continue;

                // cerchiamo START_NODE partendo da 'current_cursor' (ovvero da dove eravamo rimasti)
                auto start_it = std::find_if(
                    walk.segments.begin() + current_cursor, 
                    walk.segments.end(),
                    [&](const auto& seg) { return seg.first == start_node; }
                );

                // Se non troviamo lo start DOPO il cursore attuale, questa walk non "partecipa"
                // a questa bolla (oppure le bolle non sono ordinate).
                if (start_it == walk.segments.end()) continue;

                // print found start_it for debugging
                // std::cout << "Found start_node " << start_node << " in walk " << w_id << " at position " 
                //           << std::distance(walk.segments.begin(), start_it) << "\n";

                // cerchiamo END_NODE partendo subito dopo start_node
                auto end_it = std::find_if(
                    start_it + 1, 
                    walk.segments.end(),
                    [&](const auto& seg) { return seg.first == final_node; }
                );

                // Se troviamo anche la fine, abbiamo un match
                if (end_it != walk.segments.end()) 
                {
                    // print found end_it for debugging
                    // std::cout << "Found final_node " << final_node << " in walk " << w_id << " at position " 
                    //           << std::distance(walk.segments.begin(), end_it) << "\n";

                    // ESTRAZIONE (Start -> Inside -> End)
                    std::vector<std::string> path;
                    // Pre-allocazione: distanza tra iteratori
                    path.reserve(std::distance(start_it, end_it) + 1); 
                    
                    path.push_back(w_id);
                    
                    // Salviamo soli i nodi interni
                    for (auto it = start_it + 1; it != end_it; ++it) {
                        path.push_back(it->first);
                        // Debug print
                        // std::cout << "  Inside node: " << it->first << "\n";
                    }

                    if(path.size() < 2) {
                        // std::cout << "Warning: Extracted path has less than 2 nodes (only walk id?). Skipping.\n";
                        continue;
                    }
                    
                    haplotypes.emplace_back(std::move(path));

                    // std::cout << "--- Haplotype found for bubble " << bubble["id"] << " in walk " << w_id << "\n";

                    // AGGIORNAMENTO CURSORE
                    // Per la bolla successiva, per questa walk, cercheremo a partire dall'end_node attuale
                    current_cursor = std::distance(walk.segments.begin(), end_it);
                }
            }
            
            // continue only if haplotypes found
            if (haplotypes.empty())
                continue;

            // print haplotypes for debugging
            // std::cout << "Haplotypes found in bubble " << bubble["id"] << ":\n";
            // for (const auto &hap : haplotypes)
            // {
            //     for (size_t i = 1; i < hap.size(); ++i)
            //     {
            //         std::cout << hap[i];
            //         if (i + 1 < hap.size())
            //             std::cout << ", ";
            //     }
            //     std::cout << "\n";
            // }

            auto end_scan = std::chrono::high_resolution_clock::now();
            total_time_haplo_scan += std::chrono::duration_cast<std::chrono::microseconds>(end_scan - start_scan).count();
            // --- TIMING SECTION END ---

            struct HapKey
            {
                bool is_single;
                std::string single_node;
                std::string walk_id;
                std::vector<std::string> path_nodes;
                bool operator==(const HapKey &o) const
                {
                    return is_single == o.is_single &&
                           single_node == o.single_node &&
                           walk_id == o.walk_id &&
                           path_nodes == o.path_nodes;
                }
            };
            struct HapKeyHash
            {
                size_t operator()(const HapKey &k) const noexcept
                {
                    size_t h = hash<bool>()(k.is_single);
                    if (k.is_single)
                    {
                        h ^= hash<string>{}(k.single_node) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
                    }
                    else
                    {
                        h ^= hash<string>{}(k.walk_id) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
                        for (auto &n : k.path_nodes)
                            h ^= hash<string>{}(n) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
                    }
                    return h;
                }
            };

            // --- TIMING SECTION START: Repetition Finding (Regex & Logic) ---
            auto start_rep = std::chrono::high_resolution_clock::now();

            // Find repetitions per haplotype
            std::unordered_map<HapKey, std::vector<Repetition>, HapKeyHash> bubble_repetitions;
            for (const auto &hap : haplotypes)
            {
                if (hap.empty())
                    continue;

                if (hap.size() == 1)
                {
                    const std::string &nid = hap[0];
                    std::string sequence;
                    auto itn = gfa_data.nodes.find(nid);
                    if (itn != gfa_data.nodes.end())
                        sequence = itn->second;

                    HapKey hk;
                    hk.is_single = true;
                    hk.single_node = nid;
                    hk.walk_id = "";
                    hk.path_nodes.clear();

                    bubble_repetitions[hk] = regex_find_repetitions(sequence);
                }
                else
                {
                    std::string sequence;
                    const std::string w_id = hap[0];
                    std::vector<std::string> path_nodes;

                    auto wit = gfa_data.walks.find(w_id);
                    if (wit != gfa_data.walks.end())
                    {
                        for (const auto &seg : wit->second.segments)
                        {
                            const std::string &sid = seg.first;
                            if (std::find(hap.begin() + 1, hap.end(), sid) != hap.end())
                            {
                                path_nodes.push_back(sid);
                                auto nit = gfa_data.nodes.find(sid);
                                if (nit != gfa_data.nodes.end())
                                    sequence += nit->second;
                            }
                        }
                    }
                    else
                    {
                        for (size_t i = 1; i < hap.size(); ++i)
                        {
                            const std::string &sid = hap[i];
                            path_nodes.push_back(sid);
                            auto nit = gfa_data.nodes.find(sid);
                            if (nit != gfa_data.nodes.end())
                                sequence += nit->second;
                        }
                    }

                    HapKey hk{false, "", w_id, std::move(path_nodes)};
                    bubble_repetitions[hk] = regex_find_repetitions(sequence);
                }
            }

            // selection logic
            std::pair<Repetition, std::vector<HapKey>> selected_repetition;
            auto original_repetitions = bubble_repetitions;

            std::vector<HapKey> keys;
            keys.reserve(bubble_repetitions.size());
            for (auto &kv : bubble_repetitions)
                keys.emplace_back(kv.first);

            for (auto &node_with_rep : keys)
            {
                auto reps = bubble_repetitions[node_with_rep];
                for (auto &repetition : reps)
                {
                    const std::string &motif = get<0>(repetition);
                    if (motif.empty())
                        continue;

                    std::vector<HapKey> fusible_nodes{node_with_rep};

                    for (auto &kv2 : bubble_repetitions)
                    {
                        const auto &node = kv2.first;
                        if (node.walk_id == node_with_rep.walk_id && node.path_nodes == node_with_rep.path_nodes)
                            continue;

                        auto &vec = kv2.second;
                        vec.erase(std::remove_if(vec.begin(), vec.end(),
                                                 [&motif](const Repetition &r)
                                                 {
                                                     const auto &first = std::get<0>(r);
                                                     return first.size() > 1 && first == motif;
                                                 }),
                                  vec.end());

                        if (std::find_if(fusible_nodes.begin(), fusible_nodes.end(),
                                         [&node](const HapKey &h)
                                         { return h == node; }) == fusible_nodes.end())
                        {
                            if (node.is_single)
                            {
                                const auto &nid = node.single_node;
                                if (gfa_data.nodes.count(nid) && gfa_data.nodes[nid].find(motif) != std::string::npos)
                                {
                                    fusible_nodes.push_back(node);
                                }
                            }
                            else
                            {
                                std::string seq2;
                                for (auto &n : node.path_nodes)
                                    seq2 += gfa_data.nodes[n];
                                if (seq2.find(motif) != std::string::npos)
                                    fusible_nodes.push_back(node);
                            }
                        }
                    }

                    if (fusible_nodes.size() > 1)
                    {
                        if (selected_repetition.second.empty() ||
                            motif.size() > get<0>(selected_repetition.first).size())
                        {
                            selected_repetition = {repetition, fusible_nodes};
                        }
                    }
                }
            }

            auto end_rep = std::chrono::high_resolution_clock::now();
            total_time_rep_find += std::chrono::duration_cast<std::chrono::microseconds>(end_rep - start_rep).count();
            // --- TIMING SECTION END ---

            if (selected_repetition.second.empty())
                continue;

            // --- TIMING SECTION START: Graph Updates ---
            auto start_mod = std::chrono::high_resolution_clock::now();

            const std::string motif = get<0>(selected_repetition.first);
            unordered_map<std::string, std::string> fusible_sequences;
            for (auto &fusible : selected_repetition.second)
            {
                std::string seq;
                for (auto &n : fusible.path_nodes)
                    seq += gfa_data.nodes[n];
                fusible_sequences[fusible.walk_id] = seq;
            }

            std::string new_rep_id = new_random_id(gfa_data.nodes, "REP");
            gfa_data.nodes[new_rep_id] = motif;
            {
                Link lk;
                lk.source = new_rep_id;
                lk.source_orient = "+";
                lk.target = new_rep_id;
                lk.target_orient = "+";
                lk.overlap = "0M";
                gfa_data.links.push_back(std::move(lk));
            }

            unordered_map<std::string, std::string> up_flk_dict;
            unordered_map<std::string, std::string> dw_flk_dict;
            bool link_between_start_rep = false;
            bool link_between_rep_end = false;

            for (const auto &hk : selected_repetition.second)
            {
                const std::string w_id = hk.walk_id;
                const std::string &current_sequence = fusible_sequences[w_id];

                Repetition current_rep;
                bool found = false;
                for (auto &rep : original_repetitions[hk])
                {
                    if (std::get<0>(rep) == motif)
                    {
                        current_rep = rep;
                        found = true;
                        break;
                    }
                }
                if (!found)
                {
                    current_rep = {motif, 1, static_cast<int>(current_sequence.find(motif))};
                }

                int start_pos = std::get<2>(current_rep);
                int times = std::get<1>(current_rep);
                int end_pos = start_pos + static_cast<int>(motif.size()) * times;

                // Upstream flanking
                std::string new_up_id = "NONE";
                if (start_pos > 0)
                {
                    std::string up_seq = current_sequence.substr(0, start_pos);
                    auto itup = up_flk_dict.find(up_seq);
                    if (itup == up_flk_dict.end())
                    {
                        new_up_id = new_random_id(gfa_data.nodes, "UP_FLK");
                        up_flk_dict[up_seq] = new_up_id;
                        gfa_data.nodes[new_up_id] = up_seq;
                        {
                            Link lk1;
                            lk1.source = start_node;
                            lk1.source_orient = "+";
                            lk1.target = new_up_id;
                            lk1.target_orient = "+";
                            lk1.overlap = "0M";
                            gfa_data.links.push_back(std::move(lk1));
                        }
                        {
                            Link lk2;
                            lk2.source = new_up_id;
                            lk2.source_orient = "+";
                            lk2.target = new_rep_id;
                            lk2.target_orient = "+";
                            lk2.overlap = "0M";
                            gfa_data.links.push_back(std::move(lk2));
                        }
                    }
                    else
                    {
                        new_up_id = itup->second;
                    }
                }
                else
                {
                    if (!link_between_start_rep)
                    {
                        Link lk3;
                        lk3.source = start_node;
                        lk3.source_orient = "+";
                        lk3.target = new_rep_id;
                        lk3.target_orient = "+";
                        lk3.overlap = "0M";
                        gfa_data.links.push_back(std::move(lk3));
                        link_between_start_rep = true;
                    }
                }

                // Downstream flanking
                std::string new_dw_id = "NONE";
                if (end_pos < (int)current_sequence.size())
                {
                    std::string dw_seq = current_sequence.substr(end_pos);
                    auto itdw = dw_flk_dict.find(dw_seq);
                    if (itdw == dw_flk_dict.end())
                    {
                        new_dw_id = new_random_id(gfa_data.nodes, "DW_FLK");
                        dw_flk_dict[dw_seq] = new_dw_id;
                        gfa_data.nodes[new_dw_id] = dw_seq;
                        {
                            Link lk4;
                            lk4.source = new_rep_id;
                            lk4.source_orient = "+";
                            lk4.target = new_dw_id;
                            lk4.target_orient = "+";
                            lk4.overlap = "0M";
                            gfa_data.links.push_back(std::move(lk4));
                        }
                        {
                            Link lk5;
                            lk5.source = new_dw_id;
                            lk5.source_orient = "+";
                            lk5.target = final_node;
                            lk5.target_orient = "+";
                            lk5.overlap = "0M";
                            gfa_data.links.push_back(std::move(lk5));
                        }
                    }
                    else
                    {
                        new_dw_id = itdw->second;
                    }
                }
                else
                {
                    if (!link_between_rep_end)
                    {
                        Link lk6;
                        lk6.source = new_rep_id;
                        lk6.source_orient = "+";
                        lk6.target = final_node;
                        lk6.target_orient = "+";
                        lk6.overlap = "0M";
                        gfa_data.links.push_back(std::move(lk6));
                        link_between_rep_end = true;
                    }
                }

                // Update the corresponding WALK
                auto &segs = gfa_data.walks[w_id].segments;
                int idx = find_index_by_seg(segs, start_node);
                if (idx < 0)
                    continue;

                if (new_up_id != "NONE")
                {
                    segs.insert(segs.begin() + (idx + 1), {new_up_id, '+'});
                    idx += 1;
                }

                int insert_at = idx + 1;
                for (int i = 0; i < times; i++)
                {
                    segs.insert(segs.begin() + (insert_at + i), {new_rep_id, '+'});
                    idx += 1;
                }

                if (new_dw_id != "NONE")
                {
                    segs.insert(segs.begin() + (idx + 1), {new_dw_id, '+'});
                }
            }

            // Remove old nodes and clean walks/links
            std::vector<std::string> deleted_nodes;
            for (auto &hap : selected_repetition.second)
            {
                for (auto &n : hap.path_nodes)
                {
                    auto itn = gfa_data.nodes.find(n);
                    if (itn != gfa_data.nodes.end())
                    {
                        gfa_data.nodes.erase(n);
                        deleted_nodes.push_back(n);
                    }
                    auto &segs = gfa_data.walks[hap.walk_id].segments;
                    segs.erase(remove_if(segs.begin(), segs.end(),
                                         [&](auto &p)
                                         { return p.first == n; }),
                               segs.end());
                }
            }

            gfa_data.links.erase(remove_if(gfa_data.links.begin(), gfa_data.links.end(),
                                           [&](const Link &lk)
                                           {
                                               return find(deleted_nodes.begin(), deleted_nodes.end(), lk.source) != deleted_nodes.end() ||
                                                      find(deleted_nodes.begin(), deleted_nodes.end(), lk.target) != deleted_nodes.end();
                                           }),
                                 gfa_data.links.end());

            auto end_mod = std::chrono::high_resolution_clock::now();
            total_time_graph_mod += std::chrono::duration_cast<std::chrono::microseconds>(end_mod - start_mod).count();
            // --- TIMING SECTION END ---

        } // end for bubble
    } // end for chains

    // --- FINAL REPORT ---
    std::cerr << "--- BUBBLEGUN PROCESSING REPORT ---\n";
    std::cerr << "[info] Walk Scanning (Haplos): " << total_time_haplo_scan / 1e6 << " seconds.\n";
    std::cerr << "[info] Repetition/Regex Find:  " << total_time_rep_find / 1e6 << " seconds.\n";
    std::cerr << "[info] Graph Modification:     " << total_time_graph_mod / 1e6 << " seconds.\n";
    std::cerr << "-----------------------------------\n";
}

void write_gfa(const std::string &output_gfa, GFAData &gfa_data)
{
    std::ofstream gfa(output_gfa);
    if (!gfa)
        throw std::runtime_error("[info] Cannot open output GFA: " + output_gfa);

    // Header lines
    for (const auto &h : gfa_data.header)
        gfa << h << '\n';

    // Segments (S)
    for (const auto &[node_id, seq] : gfa_data.nodes)
        gfa << "S\t" << node_id << '\t' << seq << '\n';

    // Walks (W)
    for (const auto &[w_id, w] : gfa_data.walks)
    {
        gfa << "W\t"
            << w.sample << '\t'
            << w.hap << '\t'
            << w.seqid << '\t'
            << w.start << '\t'
            << w.end << '\t'
            << format_walk_segments(w.segments)
            << '\n';
    }

    // Links (L)
    for (const auto &lk : gfa_data.links)
        gfa << "L\t"
            << lk.source << '\t' << lk.source_orient << '\t'
            << lk.target << '\t' << lk.target_orient << '\t'
            << lk.overlap << '\n';

    gfa.close();
    if (!gfa)
        throw std::runtime_error("[info] I/O error while writing GFA: " + output_gfa);

    std::cout << "[info] New GFA written to " << output_gfa << "\n";
}

int main(int argc, char **argv)
{
    string gfa_path;
    string output_dir;
    if (!parse_args(argc, argv, gfa_path, output_dir))
    {
        usage(argv[0]);
        return 1;
    }

    // Note: the use of filesystem::path for path manipulations, which helps avoid common errors
    // with manual string concatenation and ensures platform-independent handling of file paths.
    // Also, by checking for errors immediately after attempting to create the directory, the code
    // prevents later failures when trying to write the output file.

    // Create output directory if it doesn't exist
    error_code ec;
    filesystem::create_directories(output_dir, ec);
    if (ec)
    {
        cerr << "Error: cannot create output directory '" << output_dir << "': " << ec.message() << "\n";
        return 1;
    }

    // Derive output GFA path: <output_dir>/<basename(input)>_mod.gfa
    filesystem::path gfa_p(gfa_path);
    string stem = gfa_p.stem().string(); // stem: the filename without its extension
    filesystem::path out_path = filesystem::path(output_dir) / (stem + "_mod_cpp.gfa");
    string output_gfa = out_path.string();

    // Extra variables as in the Python snippet
    string out_json = "OUT_JSON"; // name of the output JSON file that BubbleGun will create
    string fasta = "FASTA";
    string file_path = "OUT_JSON"; // TODO: why needed?
    // const int repetition_length = 1;

    cerr << "[info] Input GFA: " << gfa_path << "\n";
    cerr << "[info] Output dir: " << output_dir << "\n";
    cerr << "[info] Output GFA: " << output_gfa << "\n";

    auto start1 = std::chrono::high_resolution_clock::now();

    GFAData gfa_data = read_gfa(gfa_path); // TODO: put gfa_data with global

    auto end1 = std::chrono::high_resolution_clock::now();
    auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1);
    cerr << "[info] GFA read in " << duration1.count() / 1e6 << " seconds.\n";

    cerr << "[info] GFA loaded successfully.\n";
    cerr << "[info] Number of nodes: " << gfa_data.nodes.size() << "\n";
    cerr << "[info] Number of paths: " << gfa_data.paths.size() << "\n";
    cerr << "[info] Number of links: " << gfa_data.links.size() << "\n";
    cerr << "[info] Number of walks: " << gfa_data.walks.size() << "\n";

    string bubblegun_bin = "BubbleGun";

    auto start2 = std::chrono::high_resolution_clock::now();
    int bubblegun_rc = run_bubblegun(bubblegun_bin, gfa_path, out_json, fasta);
    auto end2 = std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2);
    cerr << "[info] BubbleGun executed in " << duration2.count() / 1e6 << " seconds.\n";

    if (bubblegun_rc != 0)
    {
        cerr << "[error] BubbleGun execution failed.\n";
        return 1;
    }

    auto start3 = std::chrono::high_resolution_clock::now();
    process_bubblegun_output(file_path, gfa_data);
    auto end3 = std::chrono::high_resolution_clock::now();
    auto duration3 = std::chrono::duration_cast<std::chrono::microseconds>(end3 - start3);
    cerr << "[info] Bubbles processed in " << duration3.count() / 1e6 << " seconds.\n";

    auto start4 = std::chrono::high_resolution_clock::now();
    write_gfa(output_gfa, gfa_data);
    auto end4 = std::chrono::high_resolution_clock::now();
    auto duration4 = std::chrono::duration_cast<std::chrono::microseconds>(end4 - start4);
    cerr << "[info] Output GFA written in " << duration4.count() / 1e6 << " seconds.\n";

    return 0;
}
