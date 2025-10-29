# bulk-RNA-seq-colab
Bulk RNA-seq analysis workflow (Salmon, tximport, edgeR/limma-voom, CAMERA/FRY) for gene expression and gene set enrichment analysis.

Open the notebook directly in Colab, or view in the GitHub project repo
----------------------------------------------------------------------------


[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/17Lot7zaajADnt8dpoGkm5evj47LnRUjL)

# Project Summary (In Silico) — Bulk RNA-seq on AVITI runs (ME vs OB)

## Abstract
We built a fully reproducible, Colab-friendly **bulk RNA-seq** pipeline that goes from **paired-end FASTQ** files to **gene-level differential expression** and **pathway enrichment**. Using four Element **AVITI** runs (two “ME” samples: SRR30769327/28; two “OB” samples: SRR30769331/32), we quantified transcripts with **Salmon** (selective-alignment), summarized to genes with **tximport**, modeled counts with **edgeR/limma-voom**, and tested **ME vs OB**. While the cohort is intentionally small (n=2 per group), coherent biology emerges: **Interferon α/γ response** is **lower in ME**, whereas **adhesion/ECM/junction**, **protein secretion**, **EMT**, and **angiogenesis** are **higher in ME**. These signals replicate across **CAMERA**, **FRY**, and **FGSEA**, and we list **leading-edge genes** that drive each pathway. The notebook saves all tables/figures and a manifest so others can **drop in their own FASTQs** and rerun end-to-end.

## Data & design
- **Instrument**: Element **AVITI**; **paired-end** bulk RNA-seq.
- **Samples**:  
  - ME: `SRR30769327`, `SRR30769328`  
  - OB: `SRR30769331`, `SRR30769332`
- **Design**: two conditions (**ME**, **OB**), 2 replicates each. No batch term in the final model (all `b1`); design can be extended later.

## Methods (concise)
1. **Quantification**  
   - **Reference**: Ensembl **GRCh38 release 111 cDNA** FASTA; built a **Salmon index** (k=31, default).  
   - **Salmon**: `--validateMappings --gcBias`; libtype autodetected (IU); one `quant.sf` per SRR under `quant/`.
2. **Gene-level import & filtering**  
   - **tximport** (length-scaled TPMs) → **edgeR** DGEList.  
   - **TMM** normalization; minimal expression filter (kept ≈13k genes).
3. **Voom + limma**  
   - `voom` precision weights; linear model `~ group`, contrast **ME − OB**; **eBayes** moderation.  
   - Multiple testing: **Benjamini–Hochberg FDR**. Optional **TREAT** for effect-size control (e.g., |log2FC|≥0.58).
4. **Diagnostics & visuals**  
   - **MDS** (voom logCPM), **volcano**, **MA**. Seaborn duplicates saved in Drive.
5. **Pathway analysis**  
   - **CAMERA** (competitive) & **FRY** (self-contained) on **MSigDB Hallmark**; barcode plots for key sets.  
   - **FGSEA (multilevel)** on a **moderated-t preranked** list (SYMBOL→ENTREZ), collections **Hallmark** and **Reactome**.  
   - **Leading-edge extraction** per significant pathway, saved as tidy CSVs.
6. **Reproducibility**  
   - `sessionInfo.txt`, `manifest.json`, consolidated `ME_vs_OB_results.xlsx`, and a zipped `results/` bundle.

## Key outputs (Drive)
- **DE tables**: `DE_ME_vs_OB.csv`, annotated `DE_ME_vs_OB_annotated.csv`, TREAT variants.  
- **Figures**: `plot_MDS_*`, `plot_volcano_*`, `plot_MA_*`, heatmaps, barcode PNGs, FGSEA curves.  
- **Pathways**: `MSigDB_Hallmark_CAMERA.csv`, `MSigDB_Hallmark_FRY.csv`, `FGSEA_Hallmark.csv`, `FGSEA_Reactome.csv`.  
- **Leading edges**: `FGSEA_leading_edge_master.csv` and per-pathway `LE_*` CSVs.

- ## Results (high-level)
- **Differential expression** (n=2 per group): modest gene-level discoveries (expected), but **coherent pathway-level signals**.
- **Hallmark FGSEA (padj < 0.05)**  
  - **Down in ME** (negative NES): *INTERFERON_ALPHA_RESPONSE*, *INTERFERON_GAMMA_RESPONSE*, “TNFα via NF-κB” (borderline).  
  - **Up in ME** (positive NES): *PROTEIN_SECRETION*, *APICAL_JUNCTION*, *EMT*, *INFLAMMATORY_RESPONSE*, *ANGIOGENESIS*.  
- **Reactome FGSEA**: **cell-cell communication**, **junction/ECM/laminin interactions**, **VEGF signaling**, **mechanotransduction** enriched in ME; **Interferon α/β signaling** depleted in ME.
- **CAMERA/FRY** agree** with FGSEA directions**, and **barcode plots** show enrichment tails in expected directions.
- **Leading-edge genes** highlight concrete drivers:  
  - Interferon-low (ME): **ISG15**, **IFITM1/3**, **OAS1**, **B2M**, **HLA-C**, **PSMB8** trending *down in ME / up in OB*.  
  - Adhesion/angiogenesis-high (ME): **CD44**, **ITGB3**, **OSMR/IL6ST**, **PLAU**, **GFPT2**, **SERPINA3**, **GPNMB** trending *up in ME*.

## Interpretation (in silico)
Together, these patterns suggest the **ME** condition is characterized by **enhanced adhesion/ECM remodeling** and **pro-angiogenic** programs with **lower basal interferon signaling** compared with **OB**. With only 2 replicates per group, we emphasize **pathway-level inference** over individual genes; nonetheless, the leading edges provide plausible markers for orthogonal validation.

## How to reuse this notebook on your own data
1. **Place files in Drive**  
   - FASTQs: `/content/drive/MyDrive/rnaseq_colab/fastq/` (paired files named `*_1.fastq.gz`, `*_2.fastq.gz`).  
   - Reference: `/content/drive/MyDrive/rnaseq_colab/ref/` (the same GRCh38.111 GTF; reuse our Salmon index or rebuild for other species).  
   - Metadata: edit `sample_sheet.csv` (columns: `sample_id,group,batch`) to match your sample folders under `quant/`.
2. **Run cells in order**  
   - **Env setup (A,B,C)** → **Index (once per reference)** → **Quant (Salmon)** → **6a/6b import** → **7a/7b voom+limma** → visuals → **pathways (CAMERA/FRY/FGSEA)**.  
   - All heavy outputs land in Drive; **runtime resets don’t lose progress**.
3. **Different species?**  
   - Download the matching **cDNA FASTA + GTF**; set `REF_DIR` accordingly; **rebuild the Salmon index**; then rerun from quant onward.
4. **More replicates or covariates**  
   - Update the design matrix (e.g., `~ group + batch`). For matched/paired designs, consider `duplicateCorrelation`.  
   - With ≥3–4 replicates per group, gene-level DE becomes more robust; pathway tests gain power.

## Limitations & safeguards
- **Small n (2 vs 2)**: individual DE genes are underpowered; prefer **pathway-level** conclusions.  
- **Batch/hidden factors** not modeled here (all `b1`); add terms if present in your metadata.  
- **Index choice** (cDNA only) is standard for transcript-level Salmon, but adding **decoys** or full genomic references can modestly affect specificity.  
- **ID mapping**: Our FGSEA uses **SYMBOL→ENTREZ** to avoid Ensembl version pitfalls.

## Reproducibility checklist
- Software versions captured in `sessionInfo.txt`.  
- Parameters and file paths recorded in `manifest.json`.  
- One-click bundle: `results/` folder + `rnaseq_results_bundle.zip`.  
- All intermediate R objects saved (`dge_ready.rds`, `voom_limma_fit.rds`) for fast resume.

## Key software
**Salmon**, **tximport**, **edgeR**, **limma-voom**, **BiocParallel**, **msigdbr**, **fgsea**, **org.Hs.eg.db**. Python visuals use **seaborn**, **matplotlib**, **scikit-learn (MDS)**.

## Conclusions (in silico)
Our Colab workflow delivers a **clean, reproducible bulk RNA-seq analysis** that starts from raw AVITI FASTQs and ends with **interpretable pathway biology** and **shareable artifacts**. Despite the minimal cohort, we observe a consistent **down-shift in interferon programs** and **up-shift in adhesion/ECM/angiogenic** axes in **ME vs OB** across **multiple, complementary** enrichment frameworks (CAMERA, FRY, FGSEA). The notebook is structured so that others can **swap in their own data**, extend the design (batch/paired), and regenerate all results and figures with minimal edits.

**Data availability & reuse**: All SRR-level outputs, DE tables, pathway summaries, leading-edge CSVs, and figures are written to Drive. Users can adjust thresholds (FDR, TREAT |log2FC|) and re-export without recomputation of alignment by reusing the saved `quant/` and RDS objects.

**Next (optional) directions**: add GO/Reactome ORA, Reactome-level enrichment visuals, report generation (HTML/PDF), and validation planning (qPCR/protein) anchored on leading-edge candidates. These are intentionally left for a future, empirical phase.


