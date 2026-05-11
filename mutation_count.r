# =============================================================
# Count synonymous and non-synonymous mutations
# from a codon-aligned FASTA file
#
# First sequence = reference
#
# Output:
# - Synonymous mutations
# - Non-synonymous mutations
# - Total mutations
# - dN/dS ratio
# - Average synonymous mutations
# - Average non-synonymous mutations
# - Average dN/dS
#
# Notes:
# - Input FASTA must be codon-aligned
# - Alignment length must preserve reading frame
# - Ambiguous bases and gaps are ignored
# =============================================================

# =============================================================
# Install packages if necessary
# =============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings", update = FALSE, ask = FALSE)
}

library(Biostrings)

# =============================================================
# Standard genetic code
# =============================================================

genetic_code <- c(

  TTT="F", TTC="F", TTA="L", TTG="L",
  TCT="S", TCC="S", TCA="S", TCG="S",
  TAT="Y", TAC="Y", TAA="*", TAG="*",
  TGT="C", TGC="C", TGA="*", TGG="W",

  CTT="L", CTC="L", CTA="L", CTG="L",
  CCT="P", CCC="P", CCA="P", CCG="P",
  CAT="H", CAC="H", CAA="Q", CAG="Q",
  CGT="R", CGC="R", CGA="R", CGG="R",

  ATT="I", ATC="I", ATA="I", ATG="M",
  ACT="T", ACC="T", ACA="T", ACG="T",
  AAT="N", AAC="N", AAA="K", AAG="K",
  AGT="S", AGC="S", AGA="R", AGG="R",

  GTT="V", GTC="V", GTA="V", GTG="V",
  GCT="A", GCC="A", GCA="A", GCG="A",
  GAT="D", GAC="D", GAA="E", GAG="E",
  GGT="G", GGC="G", GGA="G", GGG="G"
)

# =============================================================
# Codon translation function
# =============================================================

translate_codon <- function(codon) {

  aa <- genetic_code[[codon]]

  if (is.null(aa)) {
    return(NA)
  }

  return(aa)
}

# =============================================================
# Main function
# =============================================================

count_syn_nonsyn <- function(alignment_file) {

  # -----------------------------------------------------------
  # Read FASTA alignment
  # -----------------------------------------------------------

  dna <- readDNAStringSet(alignment_file)

  if (length(dna) < 2) {
    stop("FASTA file must contain at least 2 sequences.")
  }

  # -----------------------------------------------------------
  # Reference sequence
  # -----------------------------------------------------------

  ref <- dna[[1]]

  ref_len <- length(ref)

  # -----------------------------------------------------------
  # Adjust alignment length to multiple of 3
  # -----------------------------------------------------------

  remainder <- ref_len %% 3

  if (remainder != 0) {

    new_len <- ref_len - remainder

    message(
      paste(
        "Alignment length is not multiple of 3.",
        "Truncating",
        remainder,
        "base(s) from the end."
      )
    )

    dna <- subseq(dna, start = 1, end = new_len)

    ref <- dna[[1]]

    ref_len <- length(ref)
  }

  # -----------------------------------------------------------
  # Results table
  # -----------------------------------------------------------

  results <- data.frame(
    Sequence = character(),
    Synonymous = numeric(),
    NonSynonymous = numeric(),
    TotalMutations = numeric(),
    dNdS = numeric(),
    stringsAsFactors = FALSE
  )

  # ===========================================================
  # Loop through sequences
  # ===========================================================

  for (i in 2:length(dna)) {

    query <- dna[[i]]

    # ---------------------------------------------------------
    # Skip sequences with different length
    # ---------------------------------------------------------

    if (length(query) != ref_len) {

      warning(
        paste(
          "Skipping sequence:",
          names(dna)[i],
          "- different sequence length."
        )
      )

      next
    }

    syn <- 0
    nonsyn <- 0

    # =========================================================
    # Codon loop
    # =========================================================

    for (j in seq(1, ref_len, by = 3)) {

      # -------------------------------------------------------
      # Extract codons
      # -------------------------------------------------------

      ref_codon <- toupper(
        as.character(
          subseq(ref, start = j, width = 3)
        )
      )

      qry_codon <- toupper(
        as.character(
          subseq(query, start = j, width = 3)
        )
      )

      # -------------------------------------------------------
      # Ignore gaps and ambiguous bases
      # -------------------------------------------------------

      if (
        grepl("[^ACGT]", ref_codon) ||
        grepl("[^ACGT]", qry_codon)
      ) {
        next
      }

      # -------------------------------------------------------
      # Skip identical codons
      # -------------------------------------------------------

      if (ref_codon == qry_codon) {
        next
      }

      # -------------------------------------------------------
      # Translate codons
      # -------------------------------------------------------

      ref_aa <- translate_codon(ref_codon)

      qry_aa <- translate_codon(qry_codon)

      # Skip invalid codons
      if (is.na(ref_aa) || is.na(qry_aa)) {
        next
      }

      # -------------------------------------------------------
      # Count nucleotide differences
      # -------------------------------------------------------

      diffs <- sum(
        strsplit(ref_codon, "")[[1]] !=
          strsplit(qry_codon, "")[[1]]
      )

      # -------------------------------------------------------
      # Classify substitutions
      # -------------------------------------------------------

      if (ref_aa == qry_aa) {

        syn <- syn + diffs

      } else {

        nonsyn <- nonsyn + diffs
      }
    }

    # ---------------------------------------------------------
    # Total mutations
    # ---------------------------------------------------------

    total_mut <- syn + nonsyn

    # ---------------------------------------------------------
    # dN/dS calculation
    # ---------------------------------------------------------

    if (syn == 0) {

      dnds <- NA

    } else {

      dnds <- nonsyn / syn
    }

    # ---------------------------------------------------------
    # Save results
    # ---------------------------------------------------------

    results <- rbind(
      results,
      data.frame(
        Sequence = names(dna)[i],
        Synonymous = syn,
        NonSynonymous = nonsyn,
        TotalMutations = total_mut,
        dNdS = round(dnds, 4),
        stringsAsFactors = FALSE
      )
    )
  }

  # ===========================================================
  # Summary statistics
  # ===========================================================

  summary_table <- data.frame(

    Metric = c(
      "Average_Synonymous",
      "Average_NonSynonymous",
      "Average_dNdS"
    ),

    Value = c(
      round(mean(results$Synonymous, na.rm = TRUE), 4),
      round(mean(results$NonSynonymous, na.rm = TRUE), 4),
      round(mean(results$dNdS, na.rm = TRUE), 4)
    )
  )

  # ===========================================================
  # Return output
  # ===========================================================

  return(list(
    sequence_results = results,
    summary = summary_table
  ))
}

# =============================================================
# INPUT FILE
# =============================================================

# Example:
# alignment_file <- "codon_alignment.fasta"

alignment_file <- "your_alignment.fasta"

# =============================================================
# Run analysis
# =============================================================

analysis <- count_syn_nonsyn(alignment_file)

# =============================================================
# Print results
# =============================================================

cat("\n===============================\n")
cat("Sequence Results\n")
cat("===============================\n\n")

print(analysis$sequence_results)

cat("\n===============================\n")
cat("Summary Statistics\n")
cat("===============================\n\n")

print(analysis$summary)

# =============================================================
# Save output files
# =============================================================

write.csv(
  analysis$sequence_results,
  "sequence_mutation_results.csv",
  row.names = FALSE
)

write.csv(
  analysis$summary,
  "summary_mutation_statistics.csv",
  row.names = FALSE
)

cat("\nResults saved as:\n")
cat("- sequence_mutation_results.csv\n")
cat("- summary_mutation_statistics.csv\n")
