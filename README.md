# PANPHORTE - Pangenome Graphs Topology Optimizer

This repository integrates four tools to:
1. Simulate Copy Number Variations (CNV) and generate Multiple Sequence Alignments (MSA)
2. Build pangenome graphs from MSAs using VG Toolkit
3. Convert acyclic graphs to cyclic using PANPHORTE
4. Simulate reads from CNV-enriched MSAs

---

## 🔧 Prerequisites

- Python 3.8+
- [VG Toolkit](https://github.com/vgteam/vg)
- [GraphAligner](https://github.com/maickrau/GraphAligner)
- Bash (for `.sh` scripts)
- Multi-threading support (GNU `parallel` or shell threading)

---

## 🚀 Step 1: CNV_simulator (MSA Simulation)

1. **Split chromosomes**  
   ```bash
   cd CNV_simulator
   python3 fast_splitter.py
   ```  

   - Produces `length.txt` with the chromosomes lengths (put the file name inside the code).

2. **Generate CNV positions**  
   ```bash
   python3 DS_simulator_auto.py 
   ```  

   - Outputs a `.csv` per chromosome with CNV positions and copy numbers. Remember to check the path of previous file `length.txt` and the output folder.

3. **Multi-thread automation**  
   ```bash
   bash CNV_sim.sh
   ```

   - This use `CNV_sim.py` to simulate the CNV and generate the `.fasta` files per chromosome. Remember to check input and output directory and the number of threads.

---

## 🚀 Step 2: VG Toolkit (Graph Construction)

```bash
cd ..
bash vg_execution.sh
   ```

- Generates `.vg` and `.gfa` files per chromosome starting from `.fasta` files.

---

## 🚀 Step 3: PANPOHORTE (Cyclization)

1. **Convert acyclic → cyclic**  
   ```bash
   python3 main.py
   ```

2. **Multi-thread automation**  
   ```bash
   bash test2.sh
   ```

3. **Cleanup automation**  
   ```bash
   bash clear.sh
   ```

---

## 🚀 Step 4: CNV_reads_simulator (Read Simulation)

```bash
cd CNV_reads_simulator
python3 CNV_reads_sim.py
```
- Creates new haplotypes fromthe MSA with different copy numbers of the already present CNVs in `.fasta` format.

## 🚀 Step 5: Run GraphAligner to perform the alignment

```bash
cd ..
bash run_GA.sh
```
- Aligns the simulated reads `.fasta` to the graphs producing GAF files `.gaf` for each chromosome
