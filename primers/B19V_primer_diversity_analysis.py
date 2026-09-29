#!/usr/bin/env python3

"""
B19V primer/genome compatibility and diversity analysis

Inputs
------
1. Aligned B19V genomes in FASTA format
2. Degenerate primers in FASTA format

Primer names should ideally contain:
    Amplicon number + LEFT/RIGHT

Examples:
    WGS_B19V_0_LEFT
    WGS_B19V_0_RIGHT
    WGS_B19V_1_LEFT
    WGS_B19V_1_RIGHT

The script:

1. Identifies primer binding sites in a reference genome.
2. Maps the binding sites to alignment coordinates.
3. Evaluates the homologous binding site in every aligned genome.
4. Correctly handles IUPAC-degenerate primers.
5. Calculates total mismatches and mismatches in the terminal
   1, 3 and 5 nucleotides.
6. Calculates primer-level and amplicon-level compatibility.
7. Calculates genome-level panel coverage.
8. Generates CSV, Excel and PNG outputs.

IMPORTANT:
This is an IN SILICO primer compatibility analysis.
Primer mismatch does not prove experimental PCR dropout.
"""

import os
import re
import argparse
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# IUPAC DEFINITIONS
# ============================================================

IUPAC = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},

    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},

    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}

# ============================================================
# PARAMETERS
# ============================================================

TERMINAL_WINDOWS = [1, 3, 5]
MAX_ACCEPTABLE_MISMATCHES = 2
MAX_ACCEPTABLE_3P_MISMATCHES = 0

# ============================================================
# FASTA FUNCTIONS
# ============================================================

def read_fasta(path):
    sequences = {}
    current_id = None
    seq_parts = []
    
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(seq_parts).upper().replace(" ", "")
                current_id = line[1:].strip().split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
                
    if current_id is not None:
        sequences[current_id] = "".join(seq_parts).upper().replace(" ", "")
        
    if len(sequences) == 0:
        raise ValueError(f"No sequences found in {path}")
        
    return sequences

def remove_gaps(seq):
    return seq.replace("-", "").replace(".", "")

def normalize_sequence(seq):
    return re.sub(r"[^ACGTURYSWKMBDHVN\-\.]", "", seq.upper())

# ============================================================
# IUPAC MATCHING
# ============================================================

def bases_allowed_by_primer(primer_base):
    return IUPAC.get(primer_base.upper(), set())

def base_matches(primer_base, genome_base):
    primer_base = primer_base.upper()
    genome_base = genome_base.upper()

    if genome_base in {"-", ".", "N"}:
        return None

    allowed = bases_allowed_by_primer(primer_base)

    if genome_base not in {"A", "C", "G", "T"}:
        return None

    return genome_base in allowed

# ============================================================
# PRIMER ORIENTATION
# ============================================================

def reverse_complement_iupac(seq):
    complement = {
        "A": "T", "C": "G", "G": "C", "T": "A",
        "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
        "B": "V", "D": "H", "H": "D", "V": "B", "N": "N",
    }
    return "".join(complement.get(b, b) for b in seq.upper()[::-1])

def determine_orientation(primer_name):
    name = primer_name.upper()
    if "RIGHT" in name or "REV" in name or "REVERSE" in name:
        return "RIGHT"
    if "LEFT" in name or "FWD" in name or "FORWARD" in name:
        return "LEFT"
    return "UNKNOWN"

# ============================================================
# FIND PRIMER IN REFERENCE
# ============================================================

def primer_matches_sequence(primer, target):
    if len(primer) != len(target):
        raise ValueError("Primer and target must have same length")

    mismatches = []
    unknown_positions = []

    for i, (p, g) in enumerate(zip(primer, target), start=1):
        result = base_matches(p, g)
        if result is True:
            continue
        elif result is False:
            mismatches.append(i)
        else:
            unknown_positions.append(i)

    return mismatches, unknown_positions

def find_best_binding_site(primer, genome):
    genome = remove_gaps(genome)
    primer_length = len(primer)
    best = None

    for start in range(0, len(genome) - primer_length + 1):
        target = genome[start:start + primer_length]
        mismatches, unknown = primer_matches_sequence(primer, target)
        score = (len(mismatches), len(unknown))

        if best is None or score < best["score"]:
            best = {
                "start": start,
                "target": target,
                "mismatches": mismatches,
                "unknown": unknown,
                "score": score
            }
            if len(mismatches) == 0 and len(unknown) == 0:
                break

    return best

# ============================================================
# MAP UNGAPPED POSITION TO ALIGNMENT POSITION
# ============================================================

def ungapped_to_alignment_map(aligned_seq):
    mapping = []
    ungapped_position = 0
    for alignment_position, base in enumerate(aligned_seq):
        if base not in {"-", "."}:
            mapping.append((ungapped_position, alignment_position))
            ungapped_position += 1
    return mapping

def get_alignment_columns_for_binding_site(aligned_reference, ungapped_start, primer_length):
    mapping = ungapped_to_alignment_map(aligned_reference)
    if ungapped_start + primer_length > len(mapping):
        raise ValueError("Primer site exceeds reference sequence.")
    columns = [mapping[ungapped_start + i][1] for i in range(primer_length)]
    return columns

# ============================================================
# EVALUATE PRIMER AGAINST GENOME
# ============================================================

def evaluate_primer(primer_name, primer_sequence, genome_id, genome_sequence, alignment_columns):
    target = "".join(genome_sequence[col] for col in alignment_columns).upper()
    mismatches = []
    unknown = []

    for i, (p, g) in enumerate(zip(primer_sequence, target), start=1):
        result = base_matches(p, g)
        if result is False:
            mismatches.append(i)
        elif result is None:
            unknown.append(i)

    length = len(primer_sequence)
    mismatch_last_1 = [p for p in mismatches if p > length - 1]
    mismatch_last_3 = [p for p in mismatches if p > length - 3]
    mismatch_last_5 = [p for p in mismatches if p > length - 5]
    mismatch_3prime_positions = [length - p + 1 for p in mismatches]

    if len(mismatches) == 0 and len(unknown) == 0:
        status = "PERFECT"
    elif len(mismatches) == 0 and len(unknown) > 0:
        status = "UNCERTAIN"
    elif len(mismatch_last_5) > 0:
        status = "3P_RISK"
    elif len(mismatches) <= MAX_ACCEPTABLE_MISMATCHES:
        status = "LOW_RISK"
    else:
        status = "HIGH_MISMATCH"

    return {
        "Primer_Name": primer_name,
        "Genome_ID": genome_id,
        "Primer_Length": length,
        "Target_Sequence": target,
        "Total_Mismatches": len(mismatches),
        "Mismatches_last_1bp": len(mismatch_last_1),
        "Mismatches_last_3bp": len(mismatch_last_3),
        "Mismatches_last_5bp": len(mismatch_last_5),
        "Mismatch_Positions_5prime": ",".join(map(str, mismatches)),
        "Mismatch_Positions_from_3prime": ",".join(map(str, mismatch_3prime_positions)),
        "Unknown_Positions": ",".join(map(str, unknown)),
        "Unknown_Count": len(unknown),
        "Status": status
    }

# ============================================================
# PARSE AMPLICON NUMBER
# ============================================================

def get_amplicon_number(primer_name):
    patterns = [r"_(\d+)_LEFT", r"_(\d+)_RIGHT", r"_(\d+)_FWD", r"_(\d+)_REV"]
    for pattern in patterns:
        match = re.search(pattern, primer_name.upper())
        if match:
            return int(match.group(1))

    numbers = re.findall(r"\d+", primer_name)
    if numbers:
        return int(numbers[-1])

    raise ValueError(f"Could not determine amplicon number from primer name: {primer_name}")

# ============================================================
# PRIMER PAR / AMPLICON ANALYSIS
# ============================================================

def classify_amplicon(left, right):
    both_perfect = (left["Status"] == "PERFECT" and right["Status"] == "PERFECT")
    both_no_3p_mismatch = (left["Mismatches_last_5bp"] == 0 and right["Mismatches_last_5bp"] == 0)
    any_3p_mismatch = (left["Mismatches_last_5bp"] > 0 or right["Mismatches_last_5bp"] > 0)
    both_3p_mismatch = (left["Mismatches_last_5bp"] > 0 and right["Mismatches_last_5bp"] > 0)
    unknown = (left["Unknown_Count"] > 0 or right["Unknown_Count"] > 0)
    total_mismatches = (left["Total_Mismatches"] + right["Total_Mismatches"])

    if both_perfect:
        status = "PERFECT"
    elif unknown:
        status = "UNCERTAIN"
    elif not any_3p_mismatch:
        status = "NO_3P_MISMATCH"
    elif both_3p_mismatch:
        status = "BOTH_3P_RISK"
    else:
        status = "ONE_3P_RISK"

    return {
        "Amplicon": left["Amplicon"],
        "Genome_ID": left["Genome_ID"],
        "Left_Primer": left["Primer_Name"],
        "Right_Primer": right["Primer_Name"],
        "Left_Total_Mismatches": left["Total_Mismatches"],
        "Right_Total_Mismatches": right["Total_Mismatches"],
        "Total_Pair_Mismatches": total_mismatches,
        "Left_3prime5_Mismatches": left["Mismatches_last_5bp"],
        "Right_3prime5_Mismatches": right["Mismatches_last_5bp"],
        "Both_Perfect": both_perfect,
        "Both_No_3prime5_Mismatch": both_no_3p_mismatch,
        "Any_3prime5_Mismatch": any_3p_mismatch,
        "Both_3prime5_Mismatch": both_3p_mismatch,
        "Uncertain": unknown,
        "Amplicon_Status": status
    }

# ============================================================
# MAIN ANALYSIS
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="B19V primer/genome diversity analysis")
    parser.add_argument("--genomes", required=True, help="Aligned B19V genome FASTA")
    parser.add_argument("--primers", required=True, help="Degenerate primer FASTA")
    parser.add_argument("--reference", default=None, help="Reference genome ID.")
    parser.add_argument("--outdir", default="B19V_primer_analysis", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    genomes = read_fasta(args.genomes)
    primers = read_fasta(args.primers)

    print("\n==========================================")
    print("B19V PRIMER / GENOME ANALYSIS")
    print("==========================================\n")
    print(f"Genomes: {len(genomes)}")
    print(f"Primers: {len(primers)}")

    alignment_lengths = {len(seq) for seq in genomes.values()}
    if len(alignment_lengths) != 1:
        raise ValueError("All genome sequences must have the same alignment length.")

    alignment_length = list(alignment_lengths)[0]
    print(f"Alignment length: {alignment_length}")

    if args.reference:
        if args.reference not in genomes:
            raise ValueError(f"Reference {args.reference} not found.")
        reference_id = args.reference
    else:
        reference_id = list(genomes.keys())[0]

    reference_aligned = genomes[reference_id]
    reference_ungapped = remove_gaps(reference_aligned)
    print(f"Reference: {reference_id}")

    primer_info = []
    for name, sequence in primers.items():
        sequence = normalize_sequence(sequence)
        orientation = determine_orientation(name)
        if orientation == "RIGHT":
            search_sequence = reverse_complement_iupac(sequence)
        else:
            search_sequence = sequence
        amplicon = get_amplicon_number(name)
        primer_info.append({
            "Primer_Name": name,
            "Original_Sequence": sequence,
            "Search_Sequence": search_sequence,
            "Orientation": orientation,
            "Amplicon": amplicon
        })

    primer_info = sorted(primer_info, key=lambda x: (x["Amplicon"], x["Orientation"]))

    reference_sites = []
    print("\nLocating primers in reference...\n")
    for p in primer_info:
        best = find_best_binding_site(p["Search_Sequence"], reference_ungapped)
        if best is None:
            raise ValueError(f"Could not locate primer {p['Primer_Name']}")
        columns = get_alignment_columns_for_binding_site(reference_aligned, best["start"], len(p["Search_Sequence"]))
        
        p["Reference_Ungapped_Start"] = best["start"]
        p["Reference_Ungapped_End"] = best["start"] + len(p["Search_Sequence"]) - 1
        p["Reference_Target"] = best["target"]
        p["Reference_Mismatches"] = len(best["mismatches"])
        p["Alignment_Start"] = min(columns)
        p["Alignment_End"] = max(columns)
        p["Alignment_Columns"] = columns
        reference_sites.append(p)

    print("\nAnalyzing primer compatibility...\n")
    primer_results = []
    for genome_id, genome_sequence in genomes.items():
        for p in reference_sites:
            result = evaluate_primer(
                primer_name=p["Primer_Name"],
                primer_sequence=p["Search_Sequence"],
                genome_id=genome_id,
                genome_sequence=genome_sequence,
                alignment_columns=p["Alignment_Columns"]
            )
            result["Amplicon"] = p["Amplicon"]
            result["Orientation"] = p["Orientation"]
            result["Reference"] = reference_id
            result["Reference_Alignment_Start"] = p["Alignment_Start"] + 1
            result["Reference_Alignment_End"] = p["Alignment_End"] + 1
            primer_results.append(result)

    primer_df = pd.DataFrame(primer_results)

    # --------------------------------------------------------
    # Detailed Mismatch Distribution per Primer
    # --------------------------------------------------------
    mismatch_counts_list = []
    for primer_name, group in primer_df.groupby("Primer_Name"):
        row = {"Primer_Name": primer_name}
        
        # Total mismatches
        for i in range(4):
            row[f"Total_Mismatches_{i}"] = (group["Total_Mismatches"] == i).sum()
        row["Total_Mismatches_4_or_more"] = (group["Total_Mismatches"] >= 4).sum()
        
        # Mismatches in the last 5 nt (3' end)
        for i in range(3):
            row[f"Mismatches_Last_5nt_{i}"] = (group["Mismatches_last_5bp"] == i).sum()
        row["Mismatches_Last_5nt_3_or_more"] = (group["Mismatches_last_5bp"] >= 3).sum()
        
        # Mismatches in the last 3 nt (3' end)
        for i in range(3):
            row[f"Mismatches_Last_3nt_{i}"] = (group["Mismatches_last_3bp"] == i).sum()
        row["Mismatches_Last_3nt_3_or_more"] = (group["Mismatches_last_3bp"] >= 3).sum()
            
        mismatch_counts_list.append(row)
        
    primer_mismatch_counts_df = pd.DataFrame(mismatch_counts_list)


    primer_summary = []
    for primer_name, df in primer_df.groupby("Primer_Name"):
        primer_summary.append({
            "Primer_Name": primer_name,
            "Amplicon": df["Amplicon"].iloc[0],
            "Orientation": df["Orientation"].iloc[0],
            "N_Genomes": len(df),
            "Perfect_n": (df["Status"] == "PERFECT").sum(),
            "Perfect_pct": 100 * (df["Status"] == "PERFECT").mean(),
            "No_3prime5_Mismatch_n": (df["Mismatches_last_5bp"] == 0).sum(),
            "No_3prime5_Mismatch_pct": 100 * (df["Mismatches_last_5bp"] == 0).mean(),
            "Any_3prime5_Mismatch_n": (df["Mismatches_last_5bp"] > 0).sum(),
            "Any_3prime5_Mismatch_pct": 100 * (df["Mismatches_last_5bp"] > 0).mean(),
            "Any_3prime3_Mismatch_pct": 100 * (df["Mismatches_last_3bp"] > 0).mean(),
            "Any_3prime1_Mismatch_pct": 100 * (df["Mismatches_last_1bp"] > 0).mean(),
            "Mean_Total_Mismatches": df["Total_Mismatches"].mean(),
            "Max_Total_Mismatches": df["Total_Mismatches"].max(),
            "Uncertain_pct": 100 * (df["Status"] == "UNCERTAIN").mean()
        })
    primer_summary_df = pd.DataFrame(primer_summary)

    pair_results = []
    grouped = primer_df.groupby(["Amplicon", "Genome_ID"])
    for (amplicon, genome_id), group in grouped:
        left = group[group["Orientation"] == "LEFT"]
        right = group[group["Orientation"] == "RIGHT"]
        if len(left) != 1 or len(right) != 1:
            continue
        left = left.iloc[0].to_dict()
        right = right.iloc[0].to_dict()
        left["Amplicon"] = amplicon
        right["Amplicon"] = amplicon
        pair_results.append(classify_amplicon(left, right))
    amplicon_genome_df = pd.DataFrame(pair_results)

    amplicon_summary = []
    for amplicon, df in amplicon_genome_df.groupby("Amplicon"):
        amplicon_summary.append({
            "Amplicon": amplicon,
            "N_Genomes": len(df),
            "Both_Perfect_n": df["Both_Perfect"].sum(),
            "Both_Perfect_pct": 100 * df["Both_Perfect"].mean(),
            "Both_No_3prime5_Mismatch_n": df["Both_No_3prime5_Mismatch"].sum(),
            "Both_No_3prime5_Mismatch_pct": 100 * df["Both_No_3prime5_Mismatch"].mean(),
            "Any_3prime5_Mismatch_n": df["Any_3prime5_Mismatch"].sum(),
            "Any_3prime5_Mismatch_pct": 100 * df["Any_3prime5_Mismatch"].mean(),
            "Both_3prime5_Mismatch_n": df["Both_3prime5_Mismatch"].sum(),
            "Both_3prime5_Mismatch_pct": 100 * df["Both_3prime5_Mismatch"].mean(),
            "Mean_Pair_Mismatches": df["Total_Pair_Mismatches"].mean(),
            "Max_Pair_Mismatches": df["Total_Pair_Mismatches"].max(),
            "Uncertain_pct": 100 * df["Uncertain"].mean()
        })
    amplicon_summary_df = pd.DataFrame(amplicon_summary)

    genome_summary = []
    for genome_id, df in amplicon_genome_df.groupby("Genome_ID"):
        total = len(df)
        perf = df["Both_Perfect"].sum()
        no_3p = df["Both_No_3prime5_Mismatch"].sum()
        risky = df["Both_3prime5_Mismatch"].sum()
        genome_summary.append({
            "Genome_ID": genome_id,
            "Total_Amplicons": total,
            "Amplicons_Both_Perfect": perf,
            "Amplicons_Both_Perfect_pct": 100 * perf / total,
            "Amplicons_Both_No_3prime5_Mismatch": no_3p,
            "Amplicons_Both_No_3prime5_Mismatch_pct": 100 * no_3p / total,
            "Amplicons_Both_3prime5_Mismatch": risky,
            "Amplicons_Both_3prime5_Mismatch_pct": 100 * risky / total
        })
    genome_summary_df = pd.DataFrame(genome_summary)

    overall = {
        "N_Genomes": len(genomes),
        "N_Primers": len(primers),
        "N_Amplicons": len(amplicon_genome_df["Amplicon"].unique()),
        "Alignment_Length": alignment_length,
        "Perfect_Primer_Comparisons_pct": 100 * (primer_df["Status"] == "PERFECT").mean(),
        "Primer_Comparisons_No_3prime5_Mismatch_pct": 100 * (primer_df["Mismatches_last_5bp"] == 0).mean(),
        "Both_Primers_Perfect_pct": 100 * amplicon_genome_df["Both_Perfect"].mean(),
        "Both_Primers_No_3prime5_Mismatch_pct": 100 * amplicon_genome_df["Both_No_3prime5_Mismatch"].mean(),
    }
    overall_df = pd.DataFrame([overall])

    # --------------------------------------------------------
    # Save CSVs
    # --------------------------------------------------------
    primer_df.to_csv(os.path.join(args.outdir, "01_primer_genome_results.csv"), index=False)
    primer_summary_df.to_csv(os.path.join(args.outdir, "02_primer_summary.csv"), index=False)
    amplicon_genome_df.to_csv(os.path.join(args.outdir, "03_amplicon_genome_results.csv"), index=False)
    amplicon_summary_df.to_csv(os.path.join(args.outdir, "04_amplicon_summary.csv"), index=False)
    genome_summary_df.to_csv(os.path.join(args.outdir, "05_genome_panel_coverage.csv"), index=False)
    
    # Save the new mismatch counts file
    primer_mismatch_counts_df.to_csv(os.path.join(args.outdir, "09_primer_mismatch_counts.csv"), index=False)

    # --------------------------------------------------------
    # Excel workbook
    # --------------------------------------------------------
    excel_path = os.path.join(args.outdir, "B19V_primer_diversity_analysis.xlsx")
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        overall_df.to_excel(writer, sheet_name="Overall", index=False)
        primer_mismatch_counts_df.to_excel(writer, sheet_name="Mismatch_counts_per_primer", index=False)
        primer_summary_df.to_excel(writer, sheet_name="Primer_summary", index=False)
        primer_df.to_excel(writer, sheet_name="Primer_genome", index=False)
        amplicon_summary_df.to_excel(writer, sheet_name="Amplicon_summary", index=False)
        amplicon_genome_df.to_excel(writer, sheet_name="Amplicon_genome", index=False)
        genome_summary_df.to_excel(writer, sheet_name="Genome_coverage", index=False)

    print(f"\nAnalysis complete! Results saved to {args.outdir}")

if __name__ == "__main__":
    main()
