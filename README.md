# Supporting data

Gregianini et al., *Genomic epidemiology of human parvovirus B19 in Southern Brazil, 2015-2024*.

The 47 near-complete B19V genomes generated in this study are deposited in GenBank under
accession numbers PZ316512-PZ316558. The files below are the alignments, trees, metadata and
figures behind the phylogenetic, phylodynamic and statistical results reported in the paper.

## figures

Final versions of the figures in the manuscript, in PNG (300-400 dpi) and PDF.

| file | figure |
|---|---|
| `Figure_1_epidemiology.png` / `.pdf` | Figure 1: samples tested and ELISA-positive by year, map of positive cases per municipality, Ct distribution and monthly seasonality |
| `Figure_2_genotype1_phylodynamics.png` / `.pdf` | Figure 2: root-to-tip regression, time-calibrated tree of genotype 1 and detail of the Brazilian clade |
| `Supplementary_Figure_1_genotypes.png` / `.pdf` | Supplementary Figure 1: maximum-likelihood tree of the 498 genomes, showing the three genotypes |

## phylogenetics

| file | content |
|---|---|
| `alignment_498_genomes.fasta` | 47 genomes from this study, 446 curated public genomes, the reference NC_000883.2 and four genotype reference strains (LaLi, A6, D91.1, V9), aligned with MAFFT and trimmed to the coordinates of the reference |
| `ml_498_genomes.contree` | maximum-likelihood consensus tree of the alignment above, model TIM3+F+R4 selected by BIC in ModelFinder, branch support from 1,000 ultrafast bootstrap replicates |
| `alignment_genotype1_422.fasta` | genotype 1 subset used for the time-resolved analysis, with the terminal repeats masked (positions 1-615 and 5175-5596 of the reference) |
| `ml_genotype1_422.contree` | maximum-likelihood consensus tree of the genotype 1 subset, model TIM3+F+I+R3, 1,000 ultrafast bootstrap replicates |
| `timetree_genotype1_411.nwk` | time-calibrated tree, 411 tips, after excluding tips whose root-to-tip residuals exceeded four interquartile ranges |
| `node_dates_genotype1_411.json` | node dates with 95% confidence intervals and the molecular clock rate, as estimated by TreeTime |

## auspice_b19v_genotype1.json

The dated genotype 1 build in Auspice format. It can be browsed without installing anything by
dropping the file onto https://auspice.us.

## metadata

| file | content |
|---|---|
| `metadata_public_genomes_446.tsv` | accession, collection date, country and subnational location of the curated public genomes |
| `metadata_genotype1_423.tsv` | metadata of the genotype 1 subset, including the genomes from this study, whose collection dates are given as year only |

## statistics

`stats_epi.py` reproduces the tests reported in the Statistical analysis section, from the counts
in Table 1, and `stats_epi.txt` is its output.

## software

MAFFT v7.526, AliView v1.31, IQ-TREE 2 with ModelFinder and 1,000 ultrafast bootstrap replicates,
Augur v33.1.0, TreeTime v0.11.5, and Python 3.13 with SciPy 1.17.
