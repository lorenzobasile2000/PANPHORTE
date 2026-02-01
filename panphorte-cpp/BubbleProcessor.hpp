#ifndef BUBBLE_PROCESSOR_HPP
#define BUBBLE_PROCESSOR_HPP

#include "PangenomeGraph.hpp"
#include "json.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <optional>

using json = nlohmann::json;

// Haplotype Analysis Structures

/**
 * Represents the raw topology of a Superbubble found by BubbleGun.
 * Stored using efficient NodeIds instead of strings.
 */
struct BubbleContext {
    NodeId start_node;
    NodeId end_node;
    std::vector<NodeId> inside_nodes; // Nodes strictly inside the bubble

    // Metadata useful for debugging or grouping
    int chain_id;
    int bubble_id;
};

// Parsing Logic

/**
 * Parses the JSON output from BubbleGun and converts it into efficient BubbleContexts.
 * * @param json_path Path to the bubbles.json file
 * @param graph Reference to the loaded PangenomeGraph (for ID lookup)
 * @return A vector of actionable bubbles (those with >1 inside node)
 */
inline std::vector<BubbleContext> parse_bubblegun_json(const std::string& json_path, PangenomeGraph& graph) {
    std::cout << "[info] Parsing BubbleGun JSON: " << json_path << std::endl;

    std::ifstream f(json_path);
    if (!f) {
        throw std::runtime_error("Cannot open JSON file: " + json_path);
    }

    // 1. Load entire JSON into RAM
    // Note: If the JSON is massive (GBs), we might need a SAX parser, 
    // but typically BubbleGun JSONs fit in RAM.
    json data;
    try {
        f >> data;
    } catch (const json::parse_error& e) {
        throw std::runtime_error("JSON Parse Error: " + std::string(e.what()));
    }

    std::vector<BubbleContext> bubbles;

    // Reserve memory to avoid reallocations.
    // Estimate: A few thousand bubbles is common.
    bubbles.reserve(4096); 

    int total_chains = 0;
    int total_bubbles = 0;
    int skipped_bubbles = 0;

    // 2. Iterate Chains
    // Structure: { "chain_X": { "bubbles": [ ... ] }, ... }
    for (auto& [chain_key, chain_val] : data.items()) {
        total_chains++;

        if (!chain_val.contains("bubbles") || !chain_val["bubbles"].is_array()) 
            continue;

        int b_idx = 0;
        for (auto& b_obj : chain_val["bubbles"]) {
            total_bubbles++;
            b_idx++;

            // Validate structure
            if (!b_obj.contains("ends") || !b_obj.contains("inside")) {
                skipped_bubbles++;
                continue;
            }
            
            // Check array sizes
            if (b_obj["ends"].size() < 2) {
                skipped_bubbles++;
                continue;
            }

            // 3. Extract Node IDs using our Graph's interning
            try {
                std::string start_str = b_obj["ends"][0].get<std::string>();
                std::string end_str   = b_obj["ends"][1].get<std::string>();
                
                // If a node doesn't exist in our graph (e.g. filtered out during GFA load),
                // get_node_id throws. We catch and skip the bubble.
                NodeId s_id = graph.get_node_id(start_str);
                NodeId e_id = graph.get_node_id(end_str);
                
                std::vector<NodeId> inside;
                auto& in_arr = b_obj["inside"];
                inside.reserve(in_arr.size());

                bool nodes_ok = true;
                for (auto& in_node : in_arr) {
                    try {
                         inside.emplace_back(graph.get_node_id(in_node.get<std::string>()));
                    } catch (...) {
                        nodes_ok = false; 
                        break;
                    }
                }

                if (!nodes_ok) {
                    skipped_bubbles++;
                    continue;
                }
                
                // Filter: Only process bubbles with > 1 inside node (PANPHORTE logic)
                // If it is a simple bubble (1 node inside), there is no "alternative path" 
                // to optimize in the same way (or logic is trivial).
                if (inside.size() > 1) {
                    // Extract numeric chain ID from "chain_123" if possible, else 0
                    int c_id = 0;
                    try {
                         size_t underscore = chain_key.find('_');
                         if (underscore != std::string::npos) {
                             c_id = std::stoi(chain_key.substr(underscore + 1));
                         }
                    } catch(...) {}

                    bubbles.push_back({s_id, e_id, std::move(inside), c_id, b_idx});
                } else {
                    // Count as skipped or just ignore silently
                }
            } 
            catch (const std::exception& e) {
                // Usually means a node ID in the JSON was not found in the GFA
                // (e.g. BubbleGun might output nodes that were pruned or header issues)
                skipped_bubbles++;
                continue;
            }
        }
    }

    std::cout << "[info] JSON Analysis:\n"
              << "       - Chains found: " << total_chains << "\n"
              << "       - Total bubbles: " << total_bubbles << "\n"
              << "       - Skipped/Invalid: " << skipped_bubbles << "\n"
              << "       - Actionable bubbles (>1 inside nodes): " << bubbles.size() << std::endl;

    return bubbles;
}

// Analysis Logic

inline std::vector<Haplotype> extract_bubble_haplotypes(
    const BubbleContext& bubble, 
    const PangenomeGraph& graph) 
{
    std::vector<Haplotype> results;
    std::unordered_map<std::string, size_t> seq_to_hap_idx;

    // 1. Create a Set of "Allowed Inside Nodes" for O(1) lookup
    // This allows us to filter the walk efficiently.
    std::unordered_set<NodeId> inside_set;
    for (NodeId n : bubble.inside_nodes) {
        inside_set.insert(n);
    }

    // Safety check
    if (bubble.start_node >= graph.node_to_walk_indices.size()) return {};

    // 2. Iterate relevant walks
    for (size_t w_idx : graph.node_to_walk_indices[bubble.start_node]) {
        const auto& walk = graph.walks[w_idx];
        const auto& segs = walk.segments;

        // We need to find the segment of the walk that corresponds to this bubble.
        // Strategy: Find the first occurrence of Start, then look for End.
        // (Handling loops correctly requires more complex logic, but this covers 99% of cases).
        
        long start_pos = -1;
        long end_pos = -1;

        for (size_t i = 0; i < segs.size(); ++i) {
            if (segs[i].id == bubble.start_node) {
                // If we haven't found a valid start yet, mark it.
                // If we found one before but never found an end, we reset (new entry into bubble).
                start_pos = (long)i;
                end_pos = -1; // Reset end
            }
            else if (segs[i].id == bubble.end_node) {
                if (start_pos != -1) {
                    end_pos = (long)i;
                    // We found a complete Start -> ... -> End traversal.
                    // Process it immediately.
                    
                    std::string hap_seq;
                    std::vector<NodeId> hap_nodes;
                    
                    // Reserve average size to avoid reallocs
                    hap_nodes.reserve(end_pos - start_pos); 

                    // Iterate strict path between Start and End
                    for (long k = start_pos + 1; k < end_pos; ++k) {
                        NodeId nid = segs[k].id;
                        
                        hap_nodes.push_back(nid);
                        hap_seq += graph.get_sequence(nid);
                    }

                    // Store result
                    if (seq_to_hap_idx.find(hap_seq) == seq_to_hap_idx.end()) {
                        seq_to_hap_idx[hap_seq] = results.size();
                        results.push_back({hap_seq, hap_nodes, {w_idx}});
                    } else {
                        results[seq_to_hap_idx[hap_seq]].walk_indices.push_back(w_idx);
                    }

                    // Reset to search for next occurrence in same walk (if any)
                    start_pos = -1; 
                    end_pos = -1;
                    break; // Or continue if we want to handle multi-traversal
                }
            }
        }
        
        // Handle REVERSE Traversal (End -> ... -> Start)
        // BubbleGun output is Start/End agnostic relative to walk direction.
        
        long rev_start_pos = -1; 
        
        for (size_t i = 0; i < segs.size(); ++i) {
            if (segs[i].id == bubble.end_node) {
                rev_start_pos = (long)i;
            }
            else if (segs[i].id == bubble.start_node && rev_start_pos != -1) {
                long rev_end_pos = (long)i;
                
                // Found End -> Start
                std::string hap_seq;
                std::vector<NodeId> hap_nodes;

                for (long k = rev_start_pos + 1; k < rev_end_pos; ++k) {
                    NodeId nid = segs[k].id;
                    hap_nodes.push_back(nid);
                    hap_seq += graph.get_sequence(nid);
                }

                if (seq_to_hap_idx.find(hap_seq) == seq_to_hap_idx.end()) {
                    seq_to_hap_idx[hap_seq] = results.size();
                    results.push_back({hap_seq, hap_nodes, {w_idx}});
                } else {
                    results[seq_to_hap_idx[hap_seq]].walk_indices.push_back(w_idx);
                }
                
                rev_start_pos = -1; 
                break;
            }
        }
    }

    return results;
}

// Analysis Logic

inline std::optional<ModificationTask> find_common_repetition(
    const std::vector<Haplotype>& haplotypes,
    NodeId start_node,
    NodeId end_node) 
{
    if (haplotypes.size() < 2) return std::nullopt;

    std::unordered_map<size_t, std::vector<Repetition>> hap_reps;
    for (size_t i = 0; i < haplotypes.size(); ++i) {
        hap_reps[i] = PangenomeGraph::find_repetitions(haplotypes[i].sequence);
    }

    // Heuristic: Check repetitions in the first haplotype
    for (size_t i = 0; i < haplotypes.size(); ++i) {
        for (const auto& rep_a : hap_reps[i]) {
            const std::string& motif_a = std::get<0>(rep_a);
            
            int shared_count = 0;
            // Check availability in others
            for (size_t j = 0; j < haplotypes.size(); ++j) {
                if (i == j) continue;
                bool has_motif = false;
                for (const auto& rep_b : hap_reps[j]) {
                    if (std::get<0>(rep_b) == motif_a) { has_motif = true; break; }
                }
                // Also check as substring
                if (!has_motif && haplotypes[j].sequence.find(motif_a) != std::string::npos) {
                     has_motif = true;
                }
                if (has_motif) shared_count++;
            }

            if (shared_count >= 1) {
                // FOUND! Construct the task with ALL haplotypes
                ModificationTask task;
                task.bubble_start = start_node;
                task.bubble_end = end_node;
                task.haplotypes = haplotypes; // Copy all haps
                task.repeat_motif = motif_a;
                return task;
            }
        }
    }
    return std::nullopt;
}

#endif // BUBBLE_PROCESSOR_HPP