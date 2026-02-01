#ifndef PANGENOME_GRAPH_HPP
#define PANGENOME_GRAPH_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <tuple>
#include <stdexcept>
#include <string_view>
#include <charconv> // For fast number parsing if needed

using NodeId = uint32_t;

struct Segment
{
    NodeId id;
    char orientation; // '>' or '<'

    bool operator==(const Segment &other) const = default; // C++20 default comparison
};

struct Link
{
    NodeId source;
    char source_orient;
    NodeId target;
    char target_orient;
    std::string overlap;
};

struct Walk
{
    std::string sample;
    int hap;
    std::string seqid;
    std::string start;
    std::string end;

    std::vector<Segment> segments;
};

// Helper structure for the apply function
struct Haplotype
{
    std::string sequence;
    std::vector<NodeId> path_nodes;
    std::vector<size_t> walk_indices;
};

struct ModificationTask
{
    NodeId bubble_start;
    NodeId bubble_end;
    std::vector<Haplotype> haplotypes;
    std::string repeat_motif;
};

using Repetition = std::tuple<std::string, int, int>;

class PangenomeGraph
{
private:
    struct StringHash
    {
        using is_transparent = void; // Enables heterogeneous lookup

        [[nodiscard]] size_t operator()(const char *txt) const
        {
            return std::hash<std::string_view>{}(txt);
        }
        [[nodiscard]] size_t operator()(std::string_view txt) const
        {
            return std::hash<std::string_view>{}(txt);
        }
        [[nodiscard]] size_t operator()(const std::string &txt) const
        {
            return std::hash<std::string>{}(txt);
        }
    };

    NodeId id_counter = 0;

    // C++20: Uses StringHash and equal_to<> for transparent lookup
    std::unordered_map<std::string, NodeId, StringHash, std::equal_to<>> str_to_id;

    std::vector<std::string> id_to_str;
    std::vector<std::string> node_sequences;
    std::vector<bool> is_deleted;

public:
    std::vector<std::string> header;
    std::vector<Walk> walks;
    std::vector<Link> links;

    // Accepts string_view. If node exists, NO allocation happens.
    NodeId get_or_create_id(std::string_view name)
    {
        // C++20: find() accepts string_view directly
        auto it = str_to_id.find(name);
        if (it != str_to_id.end())
        {
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

    NodeId get_node_id(std::string_view name) const
    {
        auto it = str_to_id.find(name);
        if (it == str_to_id.end())
        {
            throw std::runtime_error("Node not found: " + std::string(name));
        }
        return it->second;
    }

    // Standard string interface for compatibility
    std::string get_node_name(NodeId id) const
    {
        if (id >= id_to_str.size())
            return "UNKNOWN";
        return id_to_str[id];
    }

    size_t get_node_count() const
    {
        return id_to_str.size();
    }

    void set_sequence(NodeId id, std::string_view seq)
    {
        if (id >= node_sequences.size())
        {
            node_sequences.resize(id + 1);
            is_deleted.resize(id + 1, false);
        }
        node_sequences[id] = std::string(seq);
    }

    const std::string &get_sequence(NodeId id) const
    {
        return node_sequences[id];
    }

    // only for debugging purposes
    void set_generic_id_counter(uint32_t start)
    {
        id_counter = start;
    }

    // only for debugging purposes
    std::string generate_numeric_id() {
        // Initialize on first use
        // the if won't be visited, since I call set_generic_id_counter after reading the GFA
        if (id_counter == 0) {
            // uint64_t max_val = 0;
            // for (const auto& name : id_to_str) {
            //     try {
            //         // Check if node name is purely numeric
            //         if (name.find_first_not_of("0123456789") == std::string::npos) {
            //             uint64_t val = std::stoull(name);
            //             if (val > max_val) max_val = val;
            //         }
            //     } catch(...) {}
            // }
            // id_counter = max_val + 1;
            set_generic_id_counter(static_cast<uint32_t>(get_node_count() + 1));
            std::cerr << "[debug] Initialized numeric ID counter to " << id_counter << "\n";
        } else {
            id_counter++;
        }
        return std::to_string(id_counter);
    }

    std::string generate_unique_id(std::string_view prefix)
    {
        static uint64_t counter = 0;
        std::string candidate;
        // Optimize allocation by reserving
        candidate.reserve(prefix.size() + 10);

        while (true)
        {
            counter++;
            candidate.assign(prefix);
            candidate += "_";
            candidate += std::to_string(counter);

            if (str_to_id.find(candidate) == str_to_id.end())
            {
                return candidate;
            }
        }
    }

    void load_from_gfa(const std::string &filename)
    {
        std::ifstream in(filename);
        if (!in)
            throw std::runtime_error("Cannot open GFA: " + filename);

        std::string line;
        std::vector<std::string_view> fields; // Use views for splitting

        while (std::getline(in, line))
        {
            if (line.empty() || line[0] == '#')
                continue;

            // Simple view-based split to avoid string copies during parsing
            split_view(line, fields);

            char type = line[0];
            if (type == 'S' && fields.size() >= 3)
            {
                NodeId id = get_or_create_id(fields[1]);
                set_sequence(id, fields[2]);
            }
            else if (type == 'L' && fields.size() >= 6)
            {
                Link lk;
                lk.source = get_or_create_id(fields[1]);
                lk.source_orient = (fields[2] == "+") ? '+' : '-';
                lk.target = get_or_create_id(fields[3]);
                lk.target_orient = (fields[4] == "+") ? '+' : '-';
                lk.overlap = std::string(fields[5]);
                links.push_back(std::move(lk));
            }
            else if (type == 'W' && fields.size() >= 7)
            {
                Walk w;
                w.sample = std::string(fields[1]);
                w.hap = std::stoi(std::string(fields[2])); // stoi requires string
                w.seqid = std::string(fields[3]);
                w.start = std::string(fields[4]);
                w.end = std::string(fields[5]);

                // Zero-Copy Parsing Here
                w.segments = parse_walk_line(fields[6]);

                walks.push_back(std::move(w));
            }
            else if (type == 'H')
            {
                header.push_back(line);
            }
        }
    }

    void save_to_gfa(const std::string &filename)
    {
        std::ofstream out(filename);
        if (!out)
            throw std::runtime_error("Cannot write GFA: " + filename);

        for (const auto &h : header)
            out << h << "\n";

        for (size_t i = 0; i < id_to_str.size(); ++i)
        {
            if (!is_deleted[i])
            {
                out << "S\t" << id_to_str[i] << "\t" << node_sequences[i] << "\n";
            }
        }

        for (const auto &lk : links)
        {
            if (!is_deleted[lk.source] && !is_deleted[lk.target])
            {
                out << "L\t"
                    << id_to_str[lk.source] << "\t" << lk.source_orient << "\t"
                    << id_to_str[lk.target] << "\t" << lk.target_orient << "\t"
                    << lk.overlap << "\n";
            }
        }

        for (const auto &w : walks)
        {
            out << "W\t" << w.sample << "\t" << w.hap << "\t"
                << w.seqid << "\t" << w.start << "\t" << w.end << "\t";
            out << format_walk_line(w.segments) << "\n";
        }
    }

    // Algorithms & Helpers

    // Find repetitions (Static, logic unchanged)
    static std::vector<Repetition> find_repetitions(std::string_view seq)
    {
        std::vector<Repetition> reps;
        if (seq.empty())
            return reps;

        int n = static_cast<int>(seq.size());
        const int min_rep_len = 1;

        for (int i = 0; i < n;)
        {
            bool found_rep = false;
            int max_motif_len = (n - i) / 2;

            for (int len = 1; len <= max_motif_len; ++len)
            {
                std::string_view motif = seq.substr(i, len);
                std::string_view next = seq.substr(i + len, len);

                if (motif == next)
                {
                    int count = 2;
                    int cur = i + 2 * len;
                    while (cur + len <= n && seq.substr(cur, len) == motif)
                    {
                        count++;
                        cur += len;
                    }
                    if (static_cast<int>(motif.size()) > min_rep_len)
                    {
                        reps.emplace_back(std::string(motif), count, i);
                    }
                    i = cur;
                    found_rep = true;
                    break;
                }
            }
            if (!found_rep)
                i++;
        }
        return reps;
    }

    // Replace logic (Topology)
    void replace_node_in_walks(NodeId old_id, const std::vector<Segment> &replacements)
    {
        // Since 'segments' is public, external logic can also do this,
        // but this helper ensures safety.
        for (auto &w : walks)
        {
            std::vector<Segment> new_segs;
            new_segs.reserve(w.segments.size() + 5);

            bool changed = false;
            for (const auto &step : w.segments)
            {
                if (step.id == old_id)
                {
                    // Assuming positive orientation replacement for simplicity
                    new_segs.insert(new_segs.end(), replacements.begin(), replacements.end());
                    changed = true;
                }
                else
                {
                    new_segs.push_back(step);
                }
            }
            if (changed)
                w.segments = std::move(new_segs);
        }
    }

    // Indexing

    // Maps NodeId -> List of indices in 'walks' vector (this was already done in preavious versions ideas)
    // This allows O(1) lookup of which walks visit a specific node.
    std::vector<std::vector<size_t>> node_to_walk_indices;

    void build_reverse_index()
    {
        std::cout << "[info] Building Node->Walk reverse index..." << std::endl;

        // Resize to fit all nodes
        node_to_walk_indices.assign(id_to_str.size(), {});

        // Iterate every walk and map nodes to the walk index
        for (size_t w_idx = 0; w_idx < walks.size(); ++w_idx)
        {
            const auto &walk = walks[w_idx];
            for (const auto &step : walk.segments)
            {
                // Check bounds safety
                if (step.id < node_to_walk_indices.size())
                {
                    // Optimization: Avoid duplicates if a walk loops on the same node?
                    // For simple paths, push_back is fast.
                    node_to_walk_indices[step.id].push_back(w_idx);
                }
            }
        }
        std::cout << "[info] Index built." << std::endl;
    }

    // Topology Modification (The "GTO" Step)

    void apply_bubble_optimization(const ModificationTask &task)
    {
        // 1. Create the Common Repeat Node
        // ---------------------------------------------------------------------
        // std::string rep_name = generate_unique_id("REP");
        std::string rep_name = generate_numeric_id();
        NodeId rep_id = get_or_create_id(rep_name);
        set_sequence(rep_id, task.repeat_motif);

        // Self loop for the repeat
        links.push_back({rep_id, '+', rep_id, '+', "0M"});

        // Dictionaries to deduplicate flanking nodes (Sequence -> NodeId)
        std::unordered_map<std::string, NodeId> up_cache;
        std::unordered_map<std::string, NodeId> dw_cache;

        // Track nodes to delete (Set to avoid double deletion if haps share nodes)
        std::unordered_set<NodeId> nodes_to_delete;

        // 2. Process Each Haplotype
        // ---------------------------------------------------------------------
        for (const auto &hap : task.haplotypes)
        {

            // Mark old nodes for deletion
            for (NodeId old_n : hap.path_nodes)
                nodes_to_delete.insert(old_n);

            // Analysis: Find the motif in this specific haplotype sequence
            // We re-run find_repetitions or simple string search to locate bounds
            int motif_len = (int)task.repeat_motif.size();

            // Default: Treated as one block of non-repeats if motif not found
            // (Should be rare if logic is correct, but safety first)
            int start_idx = -1;
            int count = 0;

            // Strategy: Find first occurrence of motif (Greedy)
            size_t pos = hap.sequence.find(task.repeat_motif);
            if (pos != std::string::npos)
            {
                start_idx = (int)pos;
                count = 1;
                // Check immediate repeats
                size_t next_pos = pos + motif_len;
                while (next_pos + motif_len <= hap.sequence.size())
                {
                    if (hap.sequence.substr(next_pos, motif_len) == task.repeat_motif)
                    {
                        count++;
                        next_pos += motif_len;
                    }
                    else
                    {
                        break;
                    }
                }
            }
            else
            {
                // If motif not found (e.g. deletion of the CNV), counts are 0
                start_idx = 0; // split point
                count = 0;
            }

            // 3. Define Flanks
            // -----------------------------------------------------------------
            std::string up_seq = (start_idx > 0) ? hap.sequence.substr(0, start_idx) : "";

            int end_idx = (count > 0) ? (start_idx + count * motif_len) : start_idx;
            std::string dw_seq = (end_idx < (int)hap.sequence.size()) ? hap.sequence.substr(end_idx) : "";

            // 4. Create/Get UP Flank Node
            NodeId up_id = (NodeId)-1;
            if (!up_seq.empty())
            {
                if (up_cache.find(up_seq) == up_cache.end())
                {
                    // std::string n = generate_unique_id("UP_FLK");
                    std::string n = generate_numeric_id();
                    NodeId nid = get_or_create_id(n);
                    set_sequence(nid, up_seq);
                    up_cache[up_seq] = nid;

                    // Link Start -> UP
                    links.push_back({task.bubble_start, '+', nid, '+', "0M"});
                    // Link UP -> REP
                    links.push_back({nid, '+', rep_id, '+', "0M"});
                }
                up_id = up_cache[up_seq];
            }
            else
            {
                // No UP flank: Link Start -> REP directly (if repeats exist)
                // Avoid duplicating this link if multiple empty-up haps exist
                // Ideally check existence, but GFA tolerates parallel edges.
                if (count > 0)
                {
                    links.push_back({task.bubble_start, '+', rep_id, '+', "0M"});
                }
            }

            // 5. Create/Get DW Flank Node
            NodeId dw_id = (NodeId)-1;
            if (!dw_seq.empty())
            {
                if (dw_cache.find(dw_seq) == dw_cache.end())
                {
                    // std::string n = generate_unique_id("DW_FLK");
                    std::string n = generate_numeric_id();
                    NodeId nid = get_or_create_id(n);
                    set_sequence(nid, dw_seq);
                    dw_cache[dw_seq] = nid;

                    // Link REP -> DW
                    links.push_back({rep_id, '+', nid, '+', "0M"});
                    // Link DW -> End
                    links.push_back({nid, '+', task.bubble_end, '+', "0M"});
                }
                dw_id = dw_cache[dw_seq];
            }
            else
            {
                // No DW flank: Link REP -> End directly
                if (count > 0)
                {
                    links.push_back({rep_id, '+', task.bubble_end, '+', "0M"});
                }
            }

            // Handle Case: 0 Repeats (Deletion)
            // Need link UP -> DW, or Start -> DW, or UP -> End, or Start -> End
            if (count == 0)
            {
                NodeId src = (up_id != (NodeId)-1) ? up_id : task.bubble_start;
                NodeId dst = (dw_id != (NodeId)-1) ? dw_id : task.bubble_end;
                links.push_back({src, '+', dst, '+', "0M"});
            }

            // 6. Update Walks
            // -----------------------------------------------------------------
            // Construct the replacement segment vector
            std::vector<Segment> replacement;
            if (up_id != (NodeId)-1)
                replacement.push_back({up_id, '+'});
            for (int k = 0; k < count; ++k)
                replacement.push_back({rep_id, '+'});
            if (dw_id != (NodeId)-1)
                replacement.push_back({dw_id, '+'});

            // Apply to all walks associated with this haplotype
            for (size_t w_idx : hap.walk_indices)
            {
                replace_path_in_walk(walks[w_idx], hap.path_nodes, replacement);
            }
        }

        // 7. Cleanup
        for (NodeId n : nodes_to_delete)
            delete_node(n);
    }

    /**
     * Marks a node as deleted.
     * We do not remove it from the vector to preserve the indices (IDs)
     * of all other nodes. The 'save_to_gfa' method will skip any node
     * marked as true here.
     */
    void delete_node(NodeId id)
    {
        if (id < is_deleted.size())
        {
            is_deleted[id] = true;

            // OPTIONAL: Free the memory for the sequence to save RAM immediately.
            // We keep the ID valid, but the data is gone.
            node_sequences[id].clear();
            node_sequences[id].shrink_to_fit();
        }
    }

    /**
     * Checks if a node has been deleted.
     * Useful for safety checks before modifying or traversing.
     */
    bool is_node_deleted(NodeId id) const
    {
        if (id >= is_deleted.size())
            return true; // Treat out-of-bounds as non-existent
        return is_deleted[id];
    }

private:
    /**
     * Specialized helper to replace a sequence of nodes with a new sequence.
     * Handles the fact that the path might appear multiple times or in reverse.
     */
    void replace_path_in_walk(Walk &walk, const std::vector<NodeId> &old_path, const std::vector<Segment> &new_segs)
    {
        if (old_path.empty())
            return;

        // Naive approach: Build a new segment vector
        std::vector<Segment> result;
        result.reserve(walk.segments.size() + new_segs.size());

        size_t n = walk.segments.size();
        size_t m = old_path.size();
        size_t i = 0;

        while (i < n)
        {
            // Check if old_path starts at i (Forward)
            bool match_fwd = true;
            if (i + m <= n)
            {
                for (size_t k = 0; k < m; ++k)
                {
                    if (walk.segments[i + k].id != old_path[k])
                    {
                        match_fwd = false;
                        break;
                    }
                }
            }
            else
                match_fwd = false;

            if (match_fwd)
            {
                // Insert new segments (Forward)
                result.insert(result.end(), new_segs.begin(), new_segs.end());
                i += m; // Skip old path
                continue;
            }

            // Check if old_path starts at i (Reverse)?
            // Note: BubbleGun haplotypes are extracted in fwd direction of the start node.
            // But if a walk enters the bubble in reverse (End -> Start), the IDs will be
            // the reverse of old_path.
            // For simplicity, we assume generic bubble flow.
            // Ideally, we check for reverse path here too.

            result.push_back(walk.segments[i]);
            i++;
        }
        walk.segments = std::move(result);
    }

    // Helper: Split string by tab into views (Zero-Copy)
    static void split_view(std::string_view str, std::vector<std::string_view> &out)
    {
        out.clear();
        size_t start = 0;
        size_t end = str.find('\t');
        while (end != std::string_view::npos)
        {
            out.push_back(str.substr(start, end - start));
            start = end + 1;
            end = str.find('\t', start);
        }
        out.push_back(str.substr(start));
    }

    static inline bool is_orient_char(char c) { return c == '>' || c == '<'; }

    // Zero-Copy Walk Parser using C++20 Map Lookup
    std::vector<Segment> parse_walk_line(std::string_view walk_str)
    {
        std::vector<Segment> result;
        size_t n = walk_str.size();
        result.reserve(n / 5);

        size_t i = 0;
        while (i < n && std::isspace(walk_str[i]))
            i++;

        while (i < n)
        {
            char orient = '+';
            if (walk_str[i] == '>')
            {
                orient = '+';
                i++;
            }
            else if (walk_str[i] == '<')
            {
                orient = '-';
                i++;
            }

            size_t start = i;
            while (i < n && !is_orient_char(walk_str[i]) && !std::isspace(walk_str[i]))
            {
                i++;
            }

            if (i > start)
            {
                // Creates a view directly into the line buffer
                std::string_view seg_view = walk_str.substr(start, i - start);

                // C++20 Magic: 'find' uses the view directly.
                // If ID exists, NO allocation. If new, 1 allocation.
                NodeId id = get_or_create_id(seg_view);

                result.push_back({id, orient});
            }
            while (i < n && std::isspace(walk_str[i]))
                i++;
        }
        return result;
    }

    std::string format_walk_line(const std::vector<Segment> &segs)
    {
        std::string out;
        size_t est = 0;
        for (const auto &s : segs)
            if (s.id < id_to_str.size())
                est += id_to_str[s.id].size() + 1;
        out.reserve(est);

        for (const auto &step : segs)
        {
            out.push_back(step.orientation == '+' ? '>' : '<');
            if (step.id < id_to_str.size())
                out += id_to_str[step.id];
        }
        return out;
    }
};

#endif // PANGENOME_GRAPH_HPP