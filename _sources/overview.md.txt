# === FILE: docs/overview.md ===
# Overview and Scientific Justification

## Scientific Context

The Human Immunodeficiency Virus type 1 (HIV-1) is one of the most dynamic and studied pathogens in molecular virology. Its evolution is driven by three main factors:

1.  **High Mutation Rate:** Reverse transcriptase (RT) lacks proofreading activity, generating approximately 3.4 x 10⁻⁵ substitutions per site per replication cycle.
2.  **Frequent Recombination:** The virus's ability to package two RNA genomes allows recombination when a cell is coinfected, giving rise to Circulating Recombinant Forms (CRFs) and Unique Recombinant Forms (URFs).
3.  **High Replication Rate:** Enormous daily virion production allows rapid exploration of the fitness space, facilitating immune evasion and drug resistance.

## Pipeline Workflow

The following flowchart describes the sequential stages of the HIV-1 Genomic Pipeline:

```{mermaid}
graph TD
    A[01 Download] -->|FASTA| B[02 QC]
    B -->|Filtered FASTA| C[03 Alignment]
    C -->|MSA| D[04 Model Selection]
    D -->|Best Model| E[05 Phylogeny]
    E -->|Tree File| F[06 Analysis]
    E --> G[07 Visualization]
    F -->|Reports| H[Final Results]
    G -->|Figures| H
```

## Data Flow

| Stage | Input | Output | Description |
|-------|-------|--------|-------------|
| 1. Download | Query/Config | `hiv_raw.fasta` | Retrieval of sequences from NCBI Entrez. |
| 2. QC | `hiv_raw.fasta` | `hiv_filtered.fasta`, `qc_stats.csv` | Filtering by length and ambiguity (Ns). |
| 3. Alignment | `hiv_filtered.fasta` | `hiv_aligned.fasta` | Multiple Sequence Alignment using MAFFT. |
| 4. Model Selection | `hiv_aligned.fasta` | Best Fit Model | Evolutionary model testing via IQ-TREE. |
| 5. Phylogeny | `hiv_aligned.fasta` | `hiv_phylogeny.treefile` | Maximum Likelihood tree reconstruction. |
| 6. Analysis | `hiv_filtered.fasta` | `subtypes.csv`, `resistance_report.csv` | Subtyping and drug resistance detection. |
| 7. Visualization | `hiv_phylogeny.treefile` | `hiv_tree_circular.png` | Publication-ready phylogenetic plots. |

## External Tools

| Tool | Minimum Version | Reference | Purpose |
|------|-----------------|-----------|---------|
| MAFFT | 7.4 | Katoh et al. (2013) | Large-scale MSA |
| IQ-TREE | 2.0 | Minh et al. (2020) | Phylogeny & Model selection |
| Biopython | 1.80 | Cock et al. (2009) | Sequence processing |
| Pandas | 1.5 | McKinney (2010) | Data analysis |
