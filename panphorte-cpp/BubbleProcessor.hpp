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

/**
 * Represents a pending modification to the graph.
 * This decouples the "decision" phase (Step 7) from the "mutation" phase (Step 8).
 * This structure allows Step 7 to be parallelized safely.
 */
struct ModificationTask {
    // The node that contains the repetition we want to explode
    NodeId original_node;

    // The discovered repetitive motif (e.g., "AT")
    std::string repeat_motif;

    // How many times it repeats in the specific haplotype that triggered this task
    int repeat_count;

    // Where the repetition starts (0-based index) in the original node's sequence
    int start_pos;

    // The ID of the walk/haplotype where this repetition was found
    // (Used to verify or trace back logic)
    // Optional: we might not strictly need this for the generic split logic, 
    // but it is useful for the specific logic in PANPHORTE.
    std::string source_walk_sample; 
};

// -----------------------------------------------------------------------------
// Parsing Logic
// -----------------------------------------------------------------------------

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

#endif // BUBBLE_PROCESSOR_HPP