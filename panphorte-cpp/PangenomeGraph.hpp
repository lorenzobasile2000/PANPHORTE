#ifndef PANGENOME_GRAPH_HPP
#define PANGENOME_GRAPH_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <tuple>
#include <stdexcept>
#include <string_view>
#include <charconv> // For fast number parsing if needed

// -----------------------------------------------------------------------------
// C++20 Transparent Hashing for Zero-Copy Map Lookups
// -----------------------------------------------------------------------------
struct StringHash {
    using is_transparent = void; // Enables heterogeneous lookup

    [[nodiscard]] size_t operator()(const char* txt) const {
        return std::hash<std::string_view>{}(txt);
    }
    [[nodiscard]] size_t operator()(std::string_view txt) const {
        return std::hash<std::string_view>{}(txt);
    }
    [[nodiscard]] size_t operator()(const std::string& txt) const {
        return std::hash<std::string>{}(txt);
    }
};

// -----------------------------------------------------------------------------
// Data Structures
// -----------------------------------------------------------------------------

using NodeId = uint32_t;

struct Segment {
    NodeId id;
    char orientation; // '>' or '<'

    bool operator==(const Segment& other) const = default; // C++20 default comparison
};

struct Link {
    NodeId source;
    char source_orient;
    NodeId target;
    char target_orient;
    std::string overlap;
};

struct Walk {
    std::string sample;
    int hap;
    std::string seqid;
    std::string start;
    std::string end;
    
    std::vector<Segment> segments; 
};

using Repetition = std::tuple<std::string, int, int>;

// -----------------------------------------------------------------------------
// PangenomeGraph Class
// -----------------------------------------------------------------------------

class PangenomeGraph {
private:

    // C++20: Uses StringHash and equal_to<> for transparent lookup
    std::unordered_map<std::string, NodeId, StringHash, std::equal_to<>> str_to_id;
    
    std::vector<std::string> id_to_str;
    std::vector<std::string> node_sequences;
    std::vector<bool> is_deleted;

public:
    std::vector<std::string> header;
    std::vector<Walk> walks;
    std::vector<Link> links;

    // -------------------------------------------------------------------------
    // 1. ID Management (Zero-Copy)
    // -------------------------------------------------------------------------

    // Accepts string_view. If node exists, NO allocation happens.
    NodeId get_or_create_id(std::string_view name) {
        // C++20: find() accepts string_view directly
        auto it = str_to_id.find(name);
        if (it != str_to_id.end()) {
            return it->second;
        }

        // New node: Must allocate string storage now
        NodeId new_id = static_cast<NodeId>(id_to_str.size());
        
        // Emplace constructs the string key directly from view
        str_to_id.emplace(std::string(name), new_id);
        id_to_str.emplace_back(name);
        
        node_sequences.emplace_back(""); 
        is_deleted.push_back(false);
        return new_id;
    }

    NodeId get_node_id(std::string_view name) const {
        auto it = str_to_id.find(name);
        if (it == str_to_id.end()) {
            throw std::runtime_error("Node not found: " + std::string(name));
        }
        return it->second;
    }

    // Standard string interface for compatibility
    std::string get_node_name(NodeId id) const {
        if (id >= id_to_str.size()) return "UNKNOWN";
        return id_to_str[id];
    }

    size_t get_node_count() const {
        return id_to_str.size();
    }

    void set_sequence(NodeId id, std::string_view seq) {
        if (id >= node_sequences.size()) {
            node_sequences.resize(id + 1);
            is_deleted.resize(id + 1, false);
        }
        node_sequences[id] = std::string(seq);
    }

    const std::string& get_sequence(NodeId id) const {
        return node_sequences[id];
    }

    std::string generate_unique_id(std::string_view prefix) {
        static uint64_t counter = 0;
        std::string candidate;
        // Optimize allocation by reserving
        candidate.reserve(prefix.size() + 10); 
        
        while (true) {
            counter++;
            candidate.assign(prefix);
            candidate += "_";
            candidate += std::to_string(counter);
            
            if (str_to_id.find(candidate) == str_to_id.end()) {
                return candidate;
            }
        }
    }

    // -------------------------------------------------------------------------
    // 2. I/O Methods (Optimized)
    // -------------------------------------------------------------------------

    void load_from_gfa(const std::string& filename) {
        std::ifstream in(filename);
        if (!in) throw std::runtime_error("Cannot open GFA: " + filename);

        std::string line;
        std::vector<std::string_view> fields; // Use views for splitting

        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;

            // Simple view-based split to avoid string copies during parsing
            split_view(line, fields);

            char type = line[0];
            if (type == 'S' && fields.size() >= 3) {
                NodeId id = get_or_create_id(fields[1]);
                set_sequence(id, fields[2]);
            }
            else if (type == 'L' && fields.size() >= 6) {
                Link lk;
                lk.source = get_or_create_id(fields[1]);
                lk.source_orient = (fields[2] == "+") ? '+' : '-';
                lk.target = get_or_create_id(fields[3]);
                lk.target_orient = (fields[4] == "+") ? '+' : '-';
                lk.overlap = std::string(fields[5]);
                links.push_back(std::move(lk));
            }
            else if (type == 'W' && fields.size() >= 7) {
                Walk w;
                w.sample = std::string(fields[1]);
                w.hap = std::stoi(std::string(fields[2])); // stoi requires string
                w.seqid = std::string(fields[3]);
                w.start = std::string(fields[4]);
                w.end   = std::string(fields[5]);
                
                // Zero-Copy Parsing Here
                w.segments = parse_walk_line(fields[6]);
                
                walks.push_back(std::move(w));
            }
            else if (type == 'H') {
                header.push_back(line);
            }
        }
    }

    void save_to_gfa(const std::string& filename) {
        std::ofstream out(filename);
        if (!out) throw std::runtime_error("Cannot write GFA: " + filename);

        for (const auto& h : header) out << h << "\n";

        for (size_t i = 0; i < id_to_str.size(); ++i) {
            if (!is_deleted[i]) {
                out << "S\t" << id_to_str[i] << "\t" << node_sequences[i] << "\n";
            }
        }

        for (const auto& lk : links) {
            if (!is_deleted[lk.source] && !is_deleted[lk.target]) {
                out << "L\t" 
                    << id_to_str[lk.source] << "\t" << lk.source_orient << "\t"
                    << id_to_str[lk.target] << "\t" << lk.target_orient << "\t"
                    << lk.overlap << "\n";
            }
        }

        for (const auto& w : walks) {
            out << "W\t" << w.sample << "\t" << w.hap << "\t" 
                << w.seqid << "\t" << w.start << "\t" << w.end << "\t";
            out << format_walk_line(w.segments) << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // 3. Algorithms & Helpers
    // -------------------------------------------------------------------------

    // Find repetitions (Static, logic unchanged)
    static std::vector<Repetition> find_repetitions(std::string_view seq) {
        std::vector<Repetition> reps;
        if (seq.empty()) return reps;

        int n = static_cast<int>(seq.size());
        const int min_rep_len = 1; 

        for (int i = 0; i < n; ) {
            bool found_rep = false;
            int max_motif_len = (n - i) / 2;

            for (int len = 1; len <= max_motif_len; ++len) {
                std::string_view motif = seq.substr(i, len);
                std::string_view next  = seq.substr(i + len, len);

                if (motif == next) {
                    int count = 2;
                    int cur = i + 2 * len;
                    while (cur + len <= n && seq.substr(cur, len) == motif) {
                        count++;
                        cur += len;
                    }
                    if (static_cast<int>(motif.size()) > min_rep_len) {
                         reps.emplace_back(std::string(motif), count, i);
                    }
                    i = cur; 
                    found_rep = true;
                    break; 
                }
            }
            if (!found_rep) i++;
        }
        return reps;
    }

    // Replace logic (Topology)
    void replace_node_in_walks(NodeId old_id, const std::vector<Segment>& replacements) {
        // Since 'segments' is public, external logic can also do this, 
        // but this helper ensures safety.
        for (auto& w : walks) {
            std::vector<Segment> new_segs;
            new_segs.reserve(w.segments.size() + 5); 
            
            bool changed = false;
            for (const auto& step : w.segments) {
                if (step.id == old_id) {
                    // Assuming positive orientation replacement for simplicity
                    new_segs.insert(new_segs.end(), replacements.begin(), replacements.end());
                    changed = true;
                } else {
                    new_segs.push_back(step);
                }
            }
            if (changed) w.segments = std::move(new_segs);
        }
    }

private:
    // Helper: Split string by tab into views (Zero-Copy)
    static void split_view(std::string_view str, std::vector<std::string_view>& out) {
        out.clear();
        size_t start = 0;
        size_t end = str.find('\t');
        while (end != std::string_view::npos) {
            out.push_back(str.substr(start, end - start));
            start = end + 1;
            end = str.find('\t', start);
        }
        out.push_back(str.substr(start));
    }

    static inline bool is_orient_char(char c) { return c == '>' || c == '<'; }

    // Zero-Copy Walk Parser using C++20 Map Lookup
    std::vector<Segment> parse_walk_line(std::string_view walk_str) {
        std::vector<Segment> result;
        size_t n = walk_str.size();
        result.reserve(n / 5);

        size_t i = 0;
        while (i < n && std::isspace(walk_str[i])) i++; 

        while (i < n) {
            char orient = '+';
            if (walk_str[i] == '>') { orient = '+'; i++; }
            else if (walk_str[i] == '<') { orient = '-'; i++; }

            size_t start = i;
            while (i < n && !is_orient_char(walk_str[i]) && !std::isspace(walk_str[i])) {
                i++;
            }

            if (i > start) {
                // Creates a view directly into the line buffer
                std::string_view seg_view = walk_str.substr(start, i - start);
                
                // C++20 Magic: 'find' uses the view directly.
                // If ID exists, NO allocation. If new, 1 allocation.
                NodeId id = get_or_create_id(seg_view);
                
                result.push_back({id, orient});
            }
            while (i < n && std::isspace(walk_str[i])) i++;
        }
        return result;
    }

    std::string format_walk_line(const std::vector<Segment>& segs) {
        std::string out;
        size_t est = 0;
        for(const auto& s : segs) 
            if (s.id < id_to_str.size()) est += id_to_str[s.id].size() + 1;
        out.reserve(est);

        for (const auto& step : segs) {
            out.push_back(step.orientation == '+' ? '>' : '<');
            if (step.id < id_to_str.size()) out += id_to_str[step.id];
        }
        return out;
    }
};

#endif // PANGENOME_GRAPH_HPP