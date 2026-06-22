#!/usr/bin/env python3
"""
Generate an nf-core/ampliseq samplesheet (TSV) from a directory of paired FASTQ files.

Usage:
    python make_ampliseq_samplesheet.py <fastq_dir> [--output samplesheet.tsv] [--absolute]

The sampleID is taken as the string before the first underscore in each filename.
Forward reads must contain 'R1' and reverse reads must contain 'R2' in their names.
Files are matched by their sampleID and any shared filename components.

Output columns: sampleID  forwardReads  reverseReads
"""

import argparse
import os
import re
import sys
from pathlib import Path


def find_fastq_files(directory: Path) -> list[Path]:
    """Recursively find all gzipped FASTQ files in the given directory."""
    patterns = ["*.fastq.gz", "*.fq.gz", "*.fastq", "*.fq"]
    files = []
    for pattern in patterns:
        files.extend(directory.rglob(pattern))
    return sorted(files)


def extract_sample_id(filename: str) -> str:
    """Extract the sample ID as the string before the first underscore."""
    stem = filename
    # Strip known FASTQ extensions
    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if stem.endswith(ext):
            stem = stem[: -len(ext)]
            break
    sample_id = stem.split("_")[0]
    if not sample_id:
        raise ValueError(f"Could not extract sample ID from filename: {filename}")
    return sample_id


def pair_files(files: list[Path]) -> list[tuple[str, Path, Path]]:
    """
    Pair R1 and R2 files. Returns a list of (sampleID, r1_path, r2_path) tuples.
    Matching is done by replacing 'R1' with 'R2' in the forward read filename.
    """
    r1_files = [f for f in files if re.search(r"[._]R1[._]|[._]R1$|_R1_", f.name)]
    r2_files = {f.name: f for f in files if re.search(r"[._]R2[._]|[._]R2$|_R2_", f.name)}

    pairs = []
    unmatched = []

    for r1 in r1_files:
        # Derive the expected R2 filename by substituting R1 -> R2
        r2_name = re.sub(r"([\._])R1([\._])", r"\g<1>R2\2", r1.name)
        if r2_name == r1.name:
            # Try end-of-stem substitution
            r2_name = re.sub(r"R1$", "R2", r1.stem) + "".join(r1.suffixes)

        if r2_name in r2_files:
            sample_id = extract_sample_id(r1.name)
            pairs.append((sample_id, r1, r2_files[r2_name]))
        else:
            unmatched.append(r1.name)

    if unmatched:
        print(
            f"WARNING: No R2 match found for the following R1 files:\n  "
            + "\n  ".join(unmatched),
            file=sys.stderr,
        )

    # Check for duplicate sample IDs
    ids = [p[0] for p in pairs]
    seen = set()
    duplicates = {x for x in ids if x in seen or seen.add(x)}
    if duplicates:
        print(
            f"WARNING: Duplicate sample IDs detected: {duplicates}\n"
            "Consider checking your filenames; ampliseq requires unique sampleIDs.",
            file=sys.stderr,
        )

    return sorted(pairs, key=lambda x: x[0])


def write_samplesheet(
    pairs: list[tuple[str, Path, Path]],
    output_path: Path,
    use_absolute: bool,
    fastq_dir: Path,
) -> None:
    """Write the TSV samplesheet."""
    with open(output_path, "w") as fh:
        fh.write("sampleID\tforwardReads\treverseReads\n")
        for sample_id, r1, r2 in pairs:
            if use_absolute:
                r1_str = str(r1.resolve())
                r2_str = str(r2.resolve())
            else:
                try:
                    r1_str = str(r1.relative_to(Path.cwd()))
                    r2_str = str(r2.relative_to(Path.cwd()))
                except ValueError:
                    # Fall back to absolute if relative path can't be computed
                    r1_str = str(r1.resolve())
                    r2_str = str(r2.resolve())
            fh.write(f"{sample_id}\t{r1_str}\t{r2_str}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate an nf-core/ampliseq samplesheet from a directory of paired FASTQ files."
    )
    parser.add_argument(
        "fastq_dir",
        type=Path,
        help="Directory containing paired FASTQ files (searched recursively).",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=Path("samplesheet.tsv"),
        help="Output samplesheet path (default: samplesheet.tsv).",
    )
    parser.add_argument(
        "--absolute",
        "-a",
        action="store_true",
        help="Use absolute paths in the samplesheet (default: paths relative to CWD).",
    )
    args = parser.parse_args()

    if not args.fastq_dir.is_dir():
        print(f"ERROR: '{args.fastq_dir}' is not a directory.", file=sys.stderr)
        sys.exit(1)

    files = find_fastq_files(args.fastq_dir)
    if not files:
        print(f"ERROR: No FASTQ files found in '{args.fastq_dir}'.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(files)} FASTQ file(s) in '{args.fastq_dir}'.")

    pairs = pair_files(files)
    if not pairs:
        print("ERROR: No R1/R2 pairs could be matched.", file=sys.stderr)
        sys.exit(1)

    write_samplesheet(pairs, args.output, args.absolute, args.fastq_dir)

    print(f"Wrote {len(pairs)} sample(s) to '{args.output}'.")
    print("\nFirst few rows:")
    with open(args.output) as fh:
        for i, line in enumerate(fh):
            if i > 4:
                print("  ...")
                break
            print(" ", line.rstrip())


if __name__ == "__main__":
    main()
