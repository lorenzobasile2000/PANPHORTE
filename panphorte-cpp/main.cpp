#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <filesystem>
#include <cstdio>       // popen, pclose, fgets
#include <array>
#include <cerrno>
#include <sys/wait.h>   // WIFEXITED, WEXITSTATUS (POSIX)
#include <regex>

#include <sstream>

using namespace std;

struct GFALink {
    string source;
    string source_orient;
    string target;
    string target_orient;
    string overlap;
};

struct Walk {
    std::string sample;
    int hap;
    std::string seqid;
    std::string start;
    std::string end;
    std::vector<std::pair<std::string, char>> segments; // (segId, orient) con orient in {'+','-'}
};

struct GFAData {
    vector<string> header;
    unordered_map<string, string> nodes;                 // id -> sequence
    unordered_map<string, vector<string>> paths;    // path_id -> segments (senza orientamento)
    vector<GFALink> links;

    unordered_map<string, Walk> walks;                   // walk_id -> Walk
};

// ---- Utility functions for Paths ----
static inline void split_tab(const string& s, vector<string>& out) {
    out.clear();
    stringstream ss(s);
    string item;
    while (getline(ss, item, '\t')) out.push_back(item);
}

static inline void split_comma(const string& s, vector<string>& out) {
    out.clear();
    stringstream ss(s);
    string item;
    while (getline(ss, item, ',')) out.push_back(item);
}

static inline bool ends_with_plus_or_minus(const string& s) {
    if (s.empty()) return false;
    char c = s.back();
    return c == '+' || c == '-';
}

// ---- Utility functions for Walks ----
static inline string trim(const string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

// Parse W field tipo ">s11<s12>s13" in lista di (segId, orient) con orient in {'+','-'}
static inline vector<pair<string, char>> parse_walk_string(const string& walk_str) {
    static const regex token_re(R"(([<>])([^<>\s,]+))");
    vector<pair<string, char>> result;
    smatch m;
    string s = trim(walk_str);
    auto it = s.cbegin();
    auto ed = s.cend();

    while (regex_search(it, ed, m, token_re)) {
        const string arrow = m[1].str();
        const string seg   = m[2].str();
        const char orient  = (arrow == ">") ? '+' : '-';
        result.emplace_back(seg, orient);
        it = m.suffix().first;
    }
    return result;
}

GFAData read_gfa(const string& file_gfa) {
    GFAData out;
    ifstream in(file_gfa);
    if (!in) {
        throw runtime_error("Impossibile aprire il file GFA: " + file_gfa);
    }

    string line;
    vector<string> fields;
    size_t walk_auto_id = 0; // per generare walk_id come in Python
    while (getline(in, line)) {
        if (line.empty()) continue;

        // Ignora eventuali righe di commento complete (non so se possono esserci in GFA)
        if (line[0] == '#') continue;

        char tipo_record = line[0];
        switch (tipo_record) {
            case 'H': {
                // Mantieni la riga intera come nell'implementazione Python
                out.header.push_back(line);
                break;
            }
            case 'S': {
                // Formato minimo atteso: S <sid> <sequence> ...
                split_tab(line, fields);
                if (fields.size() >= 3) {
                    const string& id_node = fields[1];
                    const string& sequence = fields[2];
                    out.nodes[id_node] = sequence;
                }
                break;
            }
            case 'P': {
                // Formato minimo atteso: P <pid> <segment(+/-),...> <overlaps or *> ...
                split_tab(line, fields);
                if (fields.size() >= 3) {
                    const string& path_id = fields[1];
                    vector<string> segs_raw;
                    split_comma(fields[2], segs_raw);

                    vector<string> segs_clean;
                    segs_clean.reserve(segs_raw.size());
                    for (const auto& s : segs_raw) {
                        if (!s.empty() && ends_with_plus_or_minus(s)) {
                            segs_clean.emplace_back(s.substr(0, s.size() - 1));
                        } else {
                            segs_clean.emplace_back(s);
                        }
                    }
                    out.paths[path_id] = std::move(segs_clean);
                }
                break;
            }
            case 'L': {
                // Formato minimo atteso: L <from> <from_orient> <to> <to_orient> <overlap> ...
                split_tab(line, fields);
                if (fields.size() >= 6) { // nel codice python è 5, ma servono almeno 6 campi
                    GFALink lk;
                    lk.source        = fields[1];
                    lk.source_orient = fields[2];
                    lk.target        = fields[3];
                    lk.target_orient = fields[4];
                    lk.overlap       = fields[5];
                    out.links.push_back(std::move(lk));
                }
                break;
            }
            case 'W': {
                // W  SampleId  HapIndex  SeqId  SeqStart  SeqEnd  Walk
                split_tab(line, fields);
                if (fields.size() >= 7) {
                    const string& sample     = fields[1];
                    const string& hap_str    = fields[2];
                    const string& seqid      = fields[3];
                    const string& seqstart   = fields[4];
                    const string& seqend     = fields[5];
                    const string& walk_field = fields[6];

                    auto segments = parse_walk_string(walk_field);

                    // walk_id = f"W_{walk_auto_id}:{sample}|{hap}|{seqid}|{seqstart}|{seqend}"
                    string walk_id;
                    walk_id.reserve(32 + sample.size() + seqid.size() + seqstart.size() + seqend.size());
                    walk_id += "W_";
                    walk_id += to_string(walk_auto_id);
                    walk_id += ":";
                    walk_id += sample;   walk_id += "|";
                    walk_id += hap_str;  walk_id += "|";
                    walk_id += seqid;    walk_id += "|";
                    walk_id += seqstart; walk_id += "|";
                    walk_id += seqend;

                    ++walk_auto_id;

                    Walk w;
                    w.sample   = sample;
                    w.hap      = stoi(hap_str);
                    w.seqid    = seqid;
                    w.start    = seqstart;
                    w.end      = seqend;
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

static void usage(const char* prog) {
    cerr << "Usage: " << prog << " -i <graph.gfa> [-o <output_dir>]\n"
            << "  -i, --input       Input GFA file (required)\n"
            << "  -o, --output_dir  Output directory for the modified GFA (default: .)\n";
}

static bool parse_args(int argc, char** argv, string& gfa_path, string& output_dir) {
    output_dir = "."; // default
    for (int i = 1; i < argc; ++i) {
        string a = argv[i];
        if ((a == "-i" || a == "--input") && i + 1 < argc) {
            gfa_path = argv[++i];
            if (gfa_path.size() < 4 || gfa_path.substr(gfa_path.size() - 4) != ".gfa") {
                cerr << "Error: input file must have .gfa extension\n";
                return false;
            }
        } else if ((a == "-o" || a == "--output_dir") && i + 1 < argc) {
        	output_dir = argv[++i];
        } else {
            cerr << "Unknown or malformed option: " << a << "\n";
            return false;
        }
    }
    return !gfa_path.empty();
}

int main(int argc, char** argv) {
    string gfa_path;
    string output_dir;
    if (!parse_args(argc, argv, gfa_path, output_dir)) {
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
    if (ec) {
        cerr << "Error: cannot create output directory '" << output_dir << "': " << ec.message() << "\n";
        return 1;
    }

    // Derive output GFA path: <output_dir>/<basename(input)>_mod.gfa
    filesystem::path gfa_p(gfa_path);
    string stem = gfa_p.stem().string(); // stem: the filename without its extension
    filesystem::path out_path = filesystem::path(output_dir) / (stem + "_mod.gfa");
    string output_gfa = out_path.string();

    // Extra variables as in the Python snippet
    string out_json = "OUT_JSON"; // name of the output JSON file that BubbleGun will create
    string fasta = "FASTA";
    string file_path = "OUT_JSON"; // TODO: why needed?
    const int repetition_length = 1;


    cerr << "[info] Input GFA: " << gfa_path << "\n";
	cerr << "[info] Output dir: " << output_dir << "\n";
    cerr << "[info] Output GFA: " << output_gfa << "\n";

    GFAData gfa_data = read_gfa(gfa_path);

    cerr << "[info] GFA loaded successfully.\n";
    cerr << "[info] Number of nodes: " << gfa_data.nodes.size() << "\n";
    cerr << "[info] Number of paths: " << gfa_data.paths.size() << "\n";
    cerr << "[info] Number of links: " << gfa_data.links.size() << "\n";

    


    return 0;
}
