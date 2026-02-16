// gfa_bubble_fusion_intids_metrics.cpp
// Compile: g++ -O3 -march=native -std=c++20 gfa_bubble_fusion_intids_metrics.cpp -o gfa_bubble_fusion -I./vendor
// Requires: nlohmann/json.hpp in include path (e.g. vendor/json.hpp)

#include <algorithm>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "json.hpp"

using json = nlohmann::json;
namespace fs = std::filesystem;

using NodeId = int32_t;

struct Segment {
  NodeId id;
  char orient; // '+' or '-'
};

struct Walk {
  std::string sample;
  int hap = 0;
  std::string seqid;
  std::string start;
  std::string end;
  std::vector<Segment> segs;
};

struct Link {
  NodeId from;
  char forient;
  NodeId to;
  char torient;
  std::string overlap;
};

struct Repeat {
  std::string motif;
  int times = 0;
  std::size_t pos = 0;
};

static inline bool parseIntSV(std::string_view s, NodeId &out) {
  if (s.empty()) return false;
  auto *b = s.data();
  auto *e = s.data() + s.size();
  NodeId val{};
  auto res = std::from_chars(b, e, val);
  if (res.ec != std::errc() || res.ptr != e) return false;
  out = val;
  return true;
}

static inline std::vector<std::string_view> splitTabs(std::string_view line) {
  std::vector<std::string_view> out;
  out.reserve(8);
  std::size_t pos = 0;
  while (pos <= line.size()) {
    std::size_t nxt = line.find('\t', pos);
    if (nxt == std::string_view::npos) {
      out.emplace_back(line.substr(pos));
      break;
    }
    out.emplace_back(line.substr(pos, nxt - pos));
    pos = nxt + 1;
  }
  return out;
}

static inline std::vector<Segment> parseWalk(std::string_view w) {
  std::vector<Segment> segs;
  segs.reserve(256);

  std::size_t i = 0;
  auto n = w.size();
  while (i < n) {
    char c = w[i];
    if (c != '>' && c != '<') { ++i; continue; }
    char orient = (c == '>') ? '+' : '-';
    std::size_t j = i + 1;
    while (j < n && w[j] != '>' && w[j] != '<') ++j;
    if (j > i + 1) {
      std::string_view tok = w.substr(i + 1, j - (i + 1));
      NodeId id{};
      if (!parseIntSV(tok, id)) {
        std::cerr << "Error: non-integer node id in W field: '" << std::string(tok) << "'\n";
        std::exit(1);
      }
      segs.push_back(Segment{id, orient});
    }
    i = j;
  }
  return segs;
}

static inline std::string formatWalk(const std::vector<Segment> &segs) {
  std::string out;
  out.reserve(segs.size() * 8);
  for (auto &s : segs) {
    out.push_back(s.orient == '+' ? '>' : '<');
    out += std::to_string(s.id);
  }
  return out;
}

// Regex-like tandem repeats: emulate re (.+?)\1+ scanning left-to-right, non-overlapping.
static inline std::vector<Repeat> findTandemRepeatsRegexLike(const std::string &s, int minMotifLen) {
  std::vector<Repeat> reps;
  const std::size_t n = s.size();
  if (n < 2) return reps;

  std::size_t pos = 0;
  while (pos + 1 < n) {
    bool found = false;

    for (std::size_t len = 1; pos + 2 * len <= n; ++len) {
      if (std::memcmp(s.data() + pos, s.data() + pos + len, len) != 0) continue;

      std::size_t times = 2;
      while (pos + (times + 1) * len <= n &&
             std::memcmp(s.data() + pos, s.data() + pos + times * len, len) == 0) {
        ++times;
      }

      if ((int)len > minMotifLen) {
        reps.push_back(Repeat{ s.substr(pos, len), (int)times, pos });
      }
      pos += times * len;
      found = true;
      break;
    }

    if (!found) ++pos;
  }

  return reps;
}

struct Haplotype {
  int walkIdx = -1;
  std::vector<NodeId> pathNodes;      // inside nodes in walk order
  std::string seq;                    // concatenated sequence of inside nodes
  std::vector<Repeat> repeats;
  std::unordered_set<std::string> repeatMotifs;
};

static inline void usage(const char *prog) {
  std::cerr
    << "Usage:\n"
    << "  " << prog << " -i <input.gfa> -o <output_dir> [--bubblegun <path>] [--repeat_min_len <int>] [--metrics <path>]\n";
}

// ---- Link dedup ----
struct LinkKey {
  NodeId from;
  char fo;
  NodeId to;
  char to_o;
  std::string ov;

  bool operator==(const LinkKey &other) const {
    return from == other.from && fo == other.fo && to == other.to && to_o == other.to_o && ov == other.ov;
  }
};

struct LinkKeyHash {
  std::size_t operator()(const LinkKey &k) const noexcept {
    std::size_t h1 = std::hash<NodeId>{}(k.from);
    std::size_t h2 = std::hash<int>{}(k.fo);
    std::size_t h3 = std::hash<NodeId>{}(k.to);
    std::size_t h4 = std::hash<int>{}(k.to_o);
    std::size_t h5 = std::hash<std::string>{}(k.ov);
    std::size_t h = h1;
    h ^= (h2 + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2));
    h ^= (h3 + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2));
    h ^= (h4 + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2));
    h ^= (h5 + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2));
    return h;
  }
};

static inline void rebuildLinkIndex(const std::vector<Link> &links,
                                    std::unordered_set<LinkKey, LinkKeyHash> &idx) {
  idx.clear();
  idx.reserve(links.size() * 2);
  for (auto &lk : links) {
    idx.insert(LinkKey{lk.from, lk.forient, lk.to, lk.torient, lk.overlap});
  }
}

static inline bool addLink(std::vector<Link> &links,
                           std::unordered_set<LinkKey, LinkKeyHash> &idx,
                           NodeId from, char fo, NodeId to, char to_o, std::string ov) {
  LinkKey key{from, fo, to, to_o, ov};
  if (idx.insert(key).second) {
    links.push_back(Link{from, fo, to, to_o, std::move(key.ov)});
    return true;
  }
  return false;
}

// ---- Next node id (incremental integer IDs) ----
static inline NodeId nextNodeId(std::unordered_map<NodeId, std::string> &nodes, NodeId &lastID) {
  while (true) {
    ++lastID;
    if (nodes.find(lastID) == nodes.end()) return lastID;
  }
}

// ---- Metrics helpers ----
struct ScalarStats {
  long long n = 0;
  long double sum = 0.0;
  long long minv = std::numeric_limits<long long>::max();
  long long maxv = std::numeric_limits<long long>::min();

  void add(long long v) {
    ++n;
    sum += (long double)v;
    minv = std::min(minv, v);
    maxv = std::max(maxv, v);
  }
  long double mean() const { return n ? (sum / (long double)n) : 0.0; }
  long long total() const { return (long long)sum; }
  bool empty() const { return n == 0; }
};

static inline std::size_t sumUniqueInsideBp(const std::vector<NodeId>& inside,
                                            const std::unordered_map<NodeId, std::string>& nodes) {
  std::size_t bp = 0;
  std::unordered_set<NodeId> seen;
  seen.reserve(inside.size() * 2);
  for (auto id : inside) {
    if (!seen.insert(id).second) continue;
    auto it = nodes.find(id);
    if (it != nodes.end()) bp += it->second.size();
  }
  return bp;
}

static inline void printStats(std::ostream& os, const std::string& name, const ScalarStats& s) {
  if (s.empty()) {
    os << name << ": (empty)\n";
    return;
  }
  os << name << ": n=" << s.n
     << " total=" << (long long)s.sum
     << " min=" << s.minv
     << " mean=" << (double)s.mean()
     << " max=" << s.maxv << "\n";
}

int main(int argc, char **argv) {
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);

  std::string inputGfa;
  std::string outputDir = ".";
  std::string bubblegun = "BubbleGun";
  int repetition_length = 1;
  std::string metricsPath; // optional override

  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if ((a == "-i" || a == "--input") && i + 1 < argc) inputGfa = argv[++i];
    else if ((a == "-o" || a == "--output_dir") && i + 1 < argc) outputDir = argv[++i];
    else if (a == "--bubblegun" && i + 1 < argc) bubblegun = argv[++i];
    else if (a == "--repeat_min_len" && i + 1 < argc) repetition_length = std::atoi(argv[++i]);
    else if (a == "--metrics" && i + 1 < argc) metricsPath = argv[++i];
    else {
      usage(argv[0]);
      return 2;
    }
  }

  if (inputGfa.empty()) {
    usage(argv[0]);
    return 2;
  }

  const fs::path inPath = fs::absolute(fs::path(inputGfa));
  const fs::path outDirAbs = fs::absolute(fs::path(outputDir));
  fs::create_directories(outDirAbs);

  const std::string base = inPath.stem().string();
  const fs::path outGfa = outDirAbs / (base + "_mod.gfa");
  const fs::path outJson = outDirAbs / (base + "_bubbles.json");
  const fs::path outFasta = outDirAbs / (base + "_bubbles.fasta");
  const fs::path outMetrics = metricsPath.empty() ? (outDirAbs / (base + "_metrics.txt")) : fs::absolute(fs::path(metricsPath));

  // ---------------- Read GFA ----------------
  std::vector<std::string> header;

  std::unordered_map<NodeId, std::string> nodes;
  nodes.reserve(1 << 20);

  std::vector<NodeId> nodeOrder;
  nodeOrder.reserve(1 << 20);

  std::vector<Walk> walks;
  walks.reserve(1 << 16);

  std::vector<Link> links;
  links.reserve(1 << 20);

  // node -> list of walks containing it (built from initial walks, used only to find candidates)
  std::unordered_map<NodeId, std::vector<int>> nodeToWalks;
  nodeToWalks.reserve(1 << 20);

  // node -> occurrences in ALL walks (kept up-to-date) for safe-delete
  std::unordered_map<NodeId, int> nodeOccInWalks;
  nodeOccInWalks.reserve(1 << 20);

  NodeId lastID = 0;

  // Initial graph totals (before modifications)
  long long initial_nodes_count = 0;
  long long initial_links_count = 0;
  long long initial_walks_count = 0;
  long long initial_graph_bp = 0;

  {
    std::ifstream in(inPath);
    if (!in) {
      std::cerr << "Error: cannot open input GFA: " << inPath << "\n";
      return 1;
    }

    std::string line;
    while (std::getline(in, line)) {
      if (line.empty()) continue;
      char t = line[0];

      if (t == 'H') {
        header.push_back(line);
      } else if (t == 'S') {
        auto cols = splitTabs(line);
        if (cols.size() >= 3) {
          NodeId id{};
          if (!parseIntSV(cols[1], id)) {
            std::cerr << "Error: S-line node id is not integer: '" << std::string(cols[1]) << "'\n";
            return 1;
          }
          nodes.emplace(id, std::string(cols[2]));
          nodeOrder.push_back(id);
          if (id > lastID) lastID = id;

          ++initial_nodes_count;
          initial_graph_bp += (long long)std::string(cols[2]).size();
        }
      } else if (t == 'L') {
        auto cols = splitTabs(line);
        if (cols.size() >= 6) {
          NodeId from{}, to{};
          if (!parseIntSV(cols[1], from) || !parseIntSV(cols[3], to)) {
            std::cerr << "Error: L-line node id is not integer\n";
            return 1;
          }
          char fo = cols[2].empty() ? '+' : cols[2][0];
          char to_o = cols[4].empty() ? '+' : cols[4][0];
          links.push_back(Link{from, fo, to, to_o, std::string(cols[5])});
          ++initial_links_count;
        }
      } else if (t == 'W') {
        auto cols = splitTabs(line);
        if (cols.size() >= 7) {
          Walk w;
          w.sample = std::string(cols[1]);
          parseIntSV(cols[2], w.hap);
          w.seqid = std::string(cols[3]);
          w.start = std::string(cols[4]);
          w.end = std::string(cols[5]);
          w.segs = parseWalk(cols[6]);

          int walkIdx = (int)walks.size();
          for (auto &s : w.segs) {
            nodeToWalks[s.id].push_back(walkIdx);
            nodeOccInWalks[s.id] += 1;
          }

          walks.push_back(std::move(w));
          ++initial_walks_count;
        }
      }
    }
  }

  std::cerr << "LastID: " << lastID << "\n";

  // Link dedup index
  std::unordered_set<LinkKey, LinkKeyHash> linkIndex;
  rebuildLinkIndex(links, linkIndex);

  // ---------------- Run BubbleGun ----------------
  {
    std::string cmd;
    cmd.reserve(512 + inPath.string().size() + outJson.string().size() + outFasta.string().size());
    cmd += bubblegun;
    cmd += " -g ";
    cmd += "\""+ inPath.string() + "\"";
    cmd += " bchains --bubble_json ";
    cmd += "\""+ outJson.string() + "\"";
    cmd += " --fasta ";
    cmd += "\""+ outFasta.string() + "\"";

    int rc = std::system(cmd.c_str());
    if (rc != 0) {
      std::cerr << "Error: BubbleGun failed with exit code: " << rc << "\n";
      return 1;
    }
  }

  // ---------------- Load bubble JSON ----------------
  json data;
  {
    std::ifstream jf(outJson);
    if (!jf) {
      std::cerr << "Error: cannot open bubble json: " << outJson << "\n";
      return 1;
    }
    jf >> data;
  }

  // ---------------- Metrics accumulators ----------------
  long long bubbles_seen = 0;
  long long bubbles_valid = 0;        // ends+inside parsed and inside>1
  long long bubbles_with_paths = 0;   // haps>=2
  long long bubbles_with_repeat = 0;  // motifs not empty
  long long bubbles_modified = 0;     // any_hap_updated==true

  ScalarStats paths_per_bubble;        // number of haplotypes considered (walk-support)
  ScalarStats unique_paths_per_bubble; // unique subpaths by sequence
  ScalarStats inside_nodes_per_bubble;
  ScalarStats inside_bp_per_bubble;

  ScalarStats motif_len_stats;         // best motif length per modified bubble
  ScalarStats support_stats;           // supporting haps per modified bubble
  ScalarStats times_stats;             // times per modified hap
  ScalarStats delta_segments_stats;    // (newSegs - oldSegs) per modified hap
  long long haps_updated_total = 0;
  long long haps_support_total = 0;
  long long start_not_found_total = 0;

  // Graph-level exact BP accounting
  long long nodes_added = 0;
  long long nodes_deleted = 0;
  long long links_added = 0;
  long long links_removed = 0;
  long long bp_nodes_added = 0;
  long long bp_nodes_deleted = 0;

  // Haplotype-level (upper bound)
  long long bp_removed_upper = 0;

  // Motif reuse
  std::unordered_map<std::string, long long> motif_use_count; // motif -> #bubbles where selected
  motif_use_count.reserve(1 << 16);

  // ---------------- Process bubbles ----------------
  if (!data.is_object()) {
    std::cerr << "Error: unexpected JSON format (expected object at root)\n";
    return 1;
  }

  for (auto &chainItem : data.items()) {
    auto &chain = chainItem.value();
    if (!chain.contains("bubbles") || !chain["bubbles"].is_array()) continue;

    for (auto &bubble : chain["bubbles"]) {
      ++bubbles_seen;

      if (!bubble.contains("ends") || !bubble["ends"].is_array() || bubble["ends"].size() < 2) continue;
      if (!bubble.contains("inside") || !bubble["inside"].is_array()) continue;

      NodeId start_node{}, final_node{};
      {
        std::string s0 = bubble["ends"][0].get<std::string>();
        std::string s1 = bubble["ends"][1].get<std::string>();
        if (!parseIntSV(s0, start_node) || !parseIntSV(s1, final_node)) {
          std::cerr << "Error: bubble ends not integer ids\n";
          return 1;
        }
      }

      // Determine order (heuristic): first walk that contains either end
      bool ok = false;
      for (auto &w : walks) {
        for (auto &s : w.segs) {
          if (s.id == start_node) { ok = true; break; }
          if (s.id == final_node) { std::swap(start_node, final_node); ok = true; break; }
        }
        if (ok) break;
      }

      std::vector<NodeId> inside;
      inside.reserve(bubble["inside"].size());
      for (auto &x : bubble["inside"]) {
        std::string sx = x.get<std::string>();
        NodeId nid{};
        if (!parseIntSV(sx, nid)) {
          std::cerr << "Error: bubble inside not integer ids\n";
          return 1;
        }
        inside.push_back(nid);
      }
      if (inside.size() <= 1) continue;

      ++bubbles_valid;
      inside_nodes_per_bubble.add((long long)inside.size());
      inside_bp_per_bubble.add((long long)sumUniqueInsideBp(inside, nodes));

      std::unordered_set<NodeId> insideSet;
      insideSet.reserve(inside.size() * 2);
      for (auto n : inside) insideSet.insert(n);

      // Candidate walks: union nodeToWalks[inside_node]
      std::unordered_set<int> candWalks;
      candWalks.reserve(inside.size() * 2);
      for (auto n : inside) {
        auto it = nodeToWalks.find(n);
        if (it == nodeToWalks.end()) continue;
        for (int wi : it->second) candWalks.insert(wi);
      }

      // Build haplotypes
      std::vector<Haplotype> haps;
      haps.reserve(candWalks.size());

      for (int wi : candWalks) {
        if (wi < 0 || wi >= (int)walks.size()) continue;
        auto &w = walks[wi];

        Haplotype h;
        h.walkIdx = wi;

        std::string seq;
        for (auto &seg : w.segs) {
          if (insideSet.find(seg.id) == insideSet.end()) continue;
          auto nit = nodes.find(seg.id);
          if (nit == nodes.end()) continue;
          h.pathNodes.push_back(seg.id);
          seq += nit->second;
        }

        if (h.pathNodes.empty()) continue;

        h.seq = std::move(seq);
        h.repeats = findTandemRepeatsRegexLike(h.seq, repetition_length);
        for (auto &r : h.repeats) h.repeatMotifs.insert(r.motif);

        haps.push_back(std::move(h));
      }

      if (haps.size() < 2) continue;

      ++bubbles_with_paths;
      paths_per_bubble.add((long long)haps.size());

      // Unique subpaths by sequence
      {
        std::unordered_set<std::string> uniq;
        uniq.reserve(haps.size() * 2);
        for (auto &h : haps) uniq.insert(h.seq);
        unique_paths_per_bubble.add((long long)uniq.size());
      }

      // Collect motif candidates
      std::vector<std::string> motifs;
      {
        std::unordered_set<std::string> seen;
        seen.reserve(256);
        for (auto &h : haps) {
          for (auto &r : h.repeats) {
            if ((int)r.motif.size() <= repetition_length) continue;
            if (seen.insert(r.motif).second) motifs.push_back(r.motif);
          }
        }
      }
      if (motifs.empty()) continue;

      ++bubbles_with_repeat;

      // Pick best motif: max length; tie-break support; then total times
      std::string bestMotif;
      int bestLen = -1;
      int bestSupport = -1;
      long long bestTimesSum = -1;
      std::vector<int> bestSupportIdx;

      for (auto &motif : motifs) {
        std::vector<int> supportIdx;
        supportIdx.reserve(haps.size());

        long long timesSum = 0;

        for (int i = 0; i < (int)haps.size(); ++i) {
          auto &h = haps[i];
          bool has = (h.repeatMotifs.find(motif) != h.repeatMotifs.end());
          if (!has) {
            if (h.seq.find(motif) != std::string::npos) has = true;
          }
          if (has) {
            supportIdx.push_back(i);
            for (auto &r : h.repeats) if (r.motif == motif) { timesSum += r.times; break; }
          }
        }

        if ((int)supportIdx.size() < 2) continue;

        int len = (int)motif.size();
        bool better =
          (len > bestLen) ||
          (len == bestLen && (int)supportIdx.size() > bestSupport) ||
          (len == bestLen && (int)supportIdx.size() == bestSupport && timesSum > bestTimesSum);

        if (better) {
          bestMotif = motif;
          bestLen = len;
          bestSupport = (int)supportIdx.size();
          bestTimesSum = timesSum;
          bestSupportIdx = std::move(supportIdx);
        }
      }

      if (bestSupportIdx.size() < 2) continue;

      // Create REP node + self-loop
      NodeId repId = nextNodeId(nodes, lastID);
      nodes.emplace(repId, bestMotif);
      nodeOrder.push_back(repId);
      nodes_added++;
      bp_nodes_added += (long long)bestMotif.size();
      if (addLink(links, linkIndex, repId, '+', repId, '+', "0M")) ++links_added;

      // motif reuse counting (per bubble where selected)
      motif_use_count[bestMotif] += 1;

      // flank dicts: sequence -> nodeId
      std::unordered_map<std::string, NodeId> upFlk;
      std::unordered_map<std::string, NodeId> dwFlk;
      upFlk.reserve(bestSupportIdx.size() * 2);
      dwFlk.reserve(bestSupportIdx.size() * 2);

      bool linkedStartRep = false;
      bool linkedRepEnd = false;

      bool any_hap_updated = false;

      // nodes we *intend* to delete (inside used by supporting haps)
      std::unordered_set<NodeId> intendedDelete;
      intendedDelete.reserve(inside.size() * 2);

      haps_support_total += (long long)bestSupportIdx.size();

      for (int hi : bestSupportIdx) {
        auto &hap = haps[hi];
        auto &w = walks[hap.walkIdx];

        // find a repeat occurrence
        std::size_t repPos = std::string::npos;
        int times = 1;

        for (auto &r : hap.repeats) {
          if (r.motif == bestMotif) {
            repPos = r.pos;
            times = r.times;
            break;
          }
        }
        if (repPos == std::string::npos) {
          repPos = hap.seq.find(bestMotif);
          times = 1;
        }
        if (repPos == std::string::npos) continue;

        // stats: times and haplotype-level upper bound
        times_stats.add((long long)times);
        if (times > 1) {
          bp_removed_upper += (long long)(times - 1) * (long long)bestMotif.size();
        }

        // Upstream flank
        std::optional<NodeId> upId;
        if (repPos > 0) {
          std::string upSeq = hap.seq.substr(0, repPos);
          auto it = upFlk.find(upSeq);
          if (it == upFlk.end()) {
            NodeId nid = nextNodeId(nodes, lastID);
            nodes.emplace(nid, upSeq);
            nodeOrder.push_back(nid);
            nodes_added++;
            bp_nodes_added += (long long)upSeq.size();
            upFlk.emplace(upSeq, nid);
            upId = nid;

            if (addLink(links, linkIndex, start_node, '+', *upId, '+', "0M")) ++links_added;
            if (addLink(links, linkIndex, *upId, '+', repId, '+', "0M")) ++links_added;
          } else {
            upId = it->second;
          }
        } else {
          if (!linkedStartRep) {
            if (addLink(links, linkIndex, start_node, '+', repId, '+', "0M")) ++links_added;
            linkedStartRep = true;
          }
        }

        // Downstream flank
        std::optional<NodeId> dwId;
        std::size_t dwPos = repPos + bestMotif.size() * (std::size_t)times;
        if (dwPos < hap.seq.size()) {
          std::string dwSeq = hap.seq.substr(dwPos);
          auto it = dwFlk.find(dwSeq);
          if (it == dwFlk.end()) {
            NodeId nid = nextNodeId(nodes, lastID);
            nodes.emplace(nid, dwSeq);
            nodeOrder.push_back(nid);
            nodes_added++;
            bp_nodes_added += (long long)dwSeq.size();
            dwFlk.emplace(dwSeq, nid);
            dwId = nid;

            if (addLink(links, linkIndex, repId, '+', *dwId, '+', "0M")) ++links_added;
            if (addLink(links, linkIndex, *dwId, '+', final_node, '+', "0M")) ++links_added;
          } else {
            dwId = it->second;
          }
        } else {
          if (!linkedRepEnd) {
            if (addLink(links, linkIndex, repId, '+', final_node, '+', "0M")) ++links_added;
            linkedRepEnd = true;
          }
        }

        // ---- Update walk (with SAFE-DELETE accounting) ----
        std::unordered_set<NodeId> delSet(hap.pathNodes.begin(), hap.pathNodes.end());

        // old size for delta
        const long long oldSegCount = (long long)w.segs.size();

        // decrement occurrence counts for old segs
        for (auto &seg : w.segs) nodeOccInWalks[seg.id] -= 1;

        std::vector<Segment> newSegs;
        newSegs.reserve(w.segs.size() + 2 + (std::size_t)times);

        bool inserted = false;
        for (auto &seg : w.segs) {
          if (delSet.find(seg.id) != delSet.end()) continue;

          newSegs.push_back(seg);

          if (!inserted && seg.id == start_node) {
            if (upId.has_value()) newSegs.push_back(Segment{*upId, '+'});
            for (int k = 0; k < times; ++k) newSegs.push_back(Segment{repId, '+'});
            if (dwId.has_value()) newSegs.push_back(Segment{*dwId, '+'});
            inserted = true;
          }
        }

        if (!inserted) {
          ++start_not_found_total;
        } else {
          any_hap_updated = true;
          ++haps_updated_total;
          delta_segments_stats.add((long long)newSegs.size() - oldSegCount);
        }

        w.segs = std::move(newSegs);

        // increment occurrence counts for new segs
        for (auto &seg : w.segs) nodeOccInWalks[seg.id] += 1;

        // mark intended deletions
        for (auto n : hap.pathNodes) intendedDelete.insert(n);
      }

      // bubble-level stats only if we actually updated at least one hap
      if (any_hap_updated) {
        ++bubbles_modified;
        motif_len_stats.add((long long)bestMotif.size());
        support_stats.add((long long)bestSupportIdx.size());
      }

      // ---- SAFE DELETE: delete only if node is NOT referenced in ANY walk (W-lines) ----
      std::unordered_set<NodeId> actuallyDeleted;
      actuallyDeleted.reserve(intendedDelete.size() * 2);

      for (auto dn : intendedDelete) {
        auto occIt = nodeOccInWalks.find(dn);
        int occ = (occIt == nodeOccInWalks.end()) ? 0 : occIt->second;
        if (occ <= 0) {
          auto itn = nodes.find(dn);
          if (itn != nodes.end()) {
            bp_nodes_deleted += (long long)itn->second.size();
            nodes_deleted++;
            nodes.erase(itn);
            actuallyDeleted.insert(dn);
          }
        }
      }

      if (!actuallyDeleted.empty()) {
        std::vector<Link> kept;
        kept.reserve(links.size());
        for (auto &lk : links) {
          bool drop = (actuallyDeleted.count(lk.from) || actuallyDeleted.count(lk.to));
          if (drop) { ++links_removed; continue; }
          kept.push_back(std::move(lk));
        }
        links = std::move(kept);
        rebuildLinkIndex(links, linkIndex);
      }
    }
  }

  // ---------------- Write output GFA ----------------
  {
    std::ofstream out(outGfa);
    if (!out) {
      std::cerr << "Error: cannot open output GFA: " << outGfa << "\n";
      return 1;
    }

    for (auto &h : header) out << h << "\n";

    for (auto id : nodeOrder) {
      auto it = nodes.find(id);
      if (it == nodes.end()) continue;
      out << "S\t" << id << "\t" << it->second << "\n";
    }

    for (auto &w : walks) {
      out << "W\t" << w.sample << "\t" << w.hap << "\t" << w.seqid
          << "\t" << w.start << "\t" << w.end << "\t" << formatWalk(w.segs) << "\n";
    }

    for (auto &lk : links) {
      out << "L\t" << lk.from << "\t" << lk.forient
          << "\t" << lk.to << "\t" << lk.torient
          << "\t" << lk.overlap << "\n";
    }
  }

  // ---------------- Write metrics (text) ----------------
  {
    // final totals
    long long final_nodes_count = 0;
    long long final_graph_bp = 0;
    for (auto &kv : nodes) {
      ++final_nodes_count;
      final_graph_bp += (long long)kv.second.size();
    }
    long long final_links_count = (long long)links.size();
    long long final_walks_count = (long long)walks.size();

    long long bp_net_saved_graph = bp_nodes_deleted - bp_nodes_added;

    std::ofstream m(outMetrics);
    if (!m) {
      std::cerr << "Error: cannot write metrics file: " << outMetrics << "\n";
      return 1;
    }

    m << "PANPHORTE / bubble fusion metrics\n";
    m << "================================\n\n";

    m << "INPUT:\n";
    m << "  input_gfa = " << inPath.string() << "\n";
    m << "  output_dir = " << outDirAbs.string() << "\n";
    m << "  bubblegun = " << bubblegun << "\n";
    m << "  repeat_min_len = " << repetition_length << "\n\n";

    m << "BUBBLES:\n";
    m << "  bubbles_seen = " << bubbles_seen << "\n";
    m << "  bubbles_valid (inside>1) = " << bubbles_valid << "\n";
    m << "  bubbles_with_paths (>=2 haplotypes) = " << bubbles_with_paths << "\n";
    m << "  bubbles_with_repeat_candidates = " << bubbles_with_repeat << "\n";
    m << "  bubbles_modified (>=1 hap updated) = " << bubbles_modified << "\n";
    if (bubbles_with_paths) {
      m << "  frac_with_repeat = " << (double)bubbles_with_repeat / (double)bubbles_with_paths << "\n";
    }
    if (bubbles_with_repeat) {
      m << "  frac_modified_given_repeat = " << (double)bubbles_modified / (double)bubbles_with_repeat << "\n";
    }
    m << "\n";

    m << "PATHS PER BUBBLE:\n";
    printStats(m, "  walk_support_paths_per_bubble", paths_per_bubble);
    printStats(m, "  unique_subpaths_per_bubble_by_seq", unique_paths_per_bubble);
    m << "\n";

    m << "BUBBLE SIZE:\n";
    printStats(m, "  inside_nodes_per_bubble", inside_nodes_per_bubble);
    printStats(m, "  inside_bp_per_bubble_unique_nodes", inside_bp_per_bubble);
    m << "\n";

    m << "REPEATS (ONLY MODIFIED BUBBLES/HAPS):\n";
    printStats(m, "  motif_len_bp", motif_len_stats);
    printStats(m, "  support_haps_per_modified_bubble", support_stats);
    printStats(m, "  repeat_times_per_modified_hap", times_stats);
    printStats(m, "  delta_segments_per_modified_hap (new-old)", delta_segments_stats);
    m << "\n";

    m << "WALK UPDATE HEALTH:\n";
    m << "  haps_support_total = " << haps_support_total << "\n";
    m << "  haps_updated_total = " << haps_updated_total << "\n";
    m << "  start_node_not_found_total = " << start_not_found_total << "\n";
    if (haps_support_total) {
      m << "  frac_start_not_found_over_support = " << (double)start_not_found_total / (double)haps_support_total << "\n";
    }
    m << "\n";

    m << "GRAPH ACCOUNTING:\n";
    m << "  initial_nodes = " << initial_nodes_count << "\n";
    m << "  final_nodes = " << final_nodes_count << "\n";
    m << "  initial_links = " << initial_links_count << "\n";
    m << "  final_links = " << final_links_count << "\n";
    m << "  initial_walks = " << initial_walks_count << "\n";
    m << "  final_walks = " << final_walks_count << "\n";
    m << "\n";

    m << "  nodes_added = " << nodes_added << "\n";
    m << "  nodes_deleted = " << nodes_deleted << "\n";
    m << "  links_added (dedup-aware) = " << links_added << "\n";
    m << "  links_removed (due to node deletes) = " << links_removed << "\n\n";

    m << "  initial_graph_bp (sum S seq lengths) = " << initial_graph_bp << "\n";
    m << "  final_graph_bp (sum S seq lengths) = " << final_graph_bp << "\n";
    m << "  bp_nodes_added = " << bp_nodes_added << "\n";
    m << "  bp_nodes_deleted = " << bp_nodes_deleted << "\n";
    m << "  bp_net_saved_graph_level (deleted - added) = " << bp_net_saved_graph << "\n";
    m << "  bp_removed_upper_haplotype_level (sum motif_len*(times-1)) = " << bp_removed_upper << "\n\n";

    m << "MOTIF REUSE:\n";
    m << "  distinct_selected_motifs = " << (long long)motif_use_count.size() << "\n";
    // top-10 motifs by usage (length only, to avoid dumping sequences huge in metrics)
    {
      std::vector<std::pair<std::string, long long>> v;
      v.reserve(motif_use_count.size());
      for (auto &kv : motif_use_count) v.push_back(kv);
      std::sort(v.begin(), v.end(), [](auto &a, auto &b){ return a.second > b.second; });

      m << "  top_motifs_by_bubble_count (count, len_bp):\n";
      for (std::size_t i = 0; i < v.size() && i < 10; ++i) {
        m << "    - " << v[i].second << " , " << v[i].first.size() << "\n";
      }
    }

    m << "\nNOTES:\n";
    m << "  - walk_support_paths_per_bubble counts how many W-walks intersect the bubble 'inside' nodes.\n";
    m << "  - unique_subpaths_per_bubble_by_seq counts unique concatenated sequences of inside nodes among those walks.\n";
    m << "  - bp_net_saved_graph_level is exact in terms of S-line sequence stored in the graph.\n";
    m << "  - bp_removed_upper_haplotype_level is an upper bound of removed redundant copies across modified haplotypes.\n";
  }

  std::cerr << "New GFA written to " << outGfa << "\n";
  std::cerr << "Metrics written to " << outMetrics << "\n";
  return 0;
}
