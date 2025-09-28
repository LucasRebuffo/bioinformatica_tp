#!/usr/bin/env python3
"""
Exercise 1 – Processing sequences from GenBank (mRNA NM_... records)

Features:
- Reads one or more GenBank files containing mRNA reference sequences (NM_...)
- For each record, generates all 6 reading frames (3 forward, 3 reverse-complement)
- Translates frames to amino acid sequences
- Detects open reading frames (ORFs) using start codon 'M' and stop '*'
- Selects the best frame as the one containing the longest ORF (configurable)
- Writes amino acid sequences to a FASTA output file, including best-frame annotation

Usage:
    python ex1.py --input NMxxxx.gbk [NMyyyy.gbk ...] --output output.fasta

If a GenBank file contains multiple records, each will be processed.

Notes:
- Requires Biopython.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable, List, Tuple

from Bio import SeqIO
from Bio.Seq import Seq


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Procesa archivos GenBank con mRNA NM_..., genera 6 marcos, "
            "traduce, detecta ORFs y escribe FASTA."
        )
    )
    parser.add_argument(
        "--input",
        "-i",
        nargs="+",
        required=True,
        help="Uno o más archivos GenBank (.gb, .gbk) con secuencias NM_...",
    )
    parser.add_argument(
        "--output",
        "-o",
        required=True,
        help="Archivo FASTA de salida con secuencias de aminoácidos",
    )
    parser.add_argument(
        "--min_orf_len",
        type=int,
        default=30,
        help=(
            "Longitud mínima (aa) para considerar un ORF. "
            "Por defecto 30 aa."
        ),
    )
    parser.add_argument(
        "--write_all_frames",
        action="store_true",
        help="Si se indica, escribe todas las traducciones de 6 marcos además de los ORFs.",
    )
    return parser.parse_args()


def generate_six_frames(dna: Seq) -> List[Tuple[str, Seq]]:
    """Return the six reading frames as (frame_label, aa_seq) without stop stripping.

    Frames are labeled: +1, +2, +3, -1, -2, -3
    """
    frames: List[Tuple[str, Seq]] = []
    # Forward frames
    for offset in range(3):
        frame_nt = dna[offset:]
        aa = frame_nt.translate(to_stop=False)
        frames.append((f"+{offset+1}", aa))
    # Reverse complement frames
    rc = dna.reverse_complement()
    for offset in range(3):
        frame_nt = rc[offset:]
        aa = frame_nt.translate(to_stop=False)
        frames.append((f"-{offset+1}", aa))
    return frames


def find_orfs_in_aa(aa_seq: Seq, min_len: int) -> List[Tuple[int, int, Seq]]:
    """Find ORFs in an amino acid sequence.

    Returns list of (start_index, end_index, orf_seq), where indices are 0-based, end exclusive.
    ORF definition: start at 'M', end at next '*' or sequence end. Must be >= min_len aa.
    """
    orfs: List[Tuple[int, int, Seq]] = []
    start_pos = None
    for i, aa in enumerate(str(aa_seq)):
        if aa == 'M' and start_pos is None:
            start_pos = i
        if aa == '*' and start_pos is not None:
            if i - start_pos >= min_len:
                orf_seq = aa_seq[start_pos:i]
                orfs.append((start_pos, i, orf_seq))
            start_pos = None
    # If sequence ends in an ORF without terminal stop, consider it as well
    if start_pos is not None:
        if len(aa_seq) - start_pos >= min_len:
            orf_seq = aa_seq[start_pos:]
            orfs.append((start_pos, len(aa_seq), orf_seq))
    return orfs


def select_best_frame(frames_orfs: List[Tuple[str, List[Tuple[int, int, Seq]]]]) -> Tuple[str, Tuple[int, int, Seq]] | None:
    """Select the best frame by longest ORF length.

    Returns (frame_label, (start, end, orf_seq)) or None if no ORFs found.
    """
    best: Tuple[str, Tuple[int, int, Seq]] | None = None
    best_len = -1
    for frame_label, orfs in frames_orfs:
        for start, end, orf_seq in orfs:
            if len(orf_seq) > best_len:
                best_len = len(orf_seq)
                best = (frame_label, (start, end, orf_seq))
    return best


def fasta_header(record_id: str, frame_label: str, start: int, end: int, is_best: bool) -> str:
    best_tag = " best_frame=true" if is_best else ""
    return f">{record_id} frame={frame_label} aa_start={start} aa_end={end}{best_tag}"


def safe_record_id(rec) -> str:
    # Prefer accession or id
    rid = getattr(rec, "id", None) or getattr(rec, "name", None) or "unknown"
    # Remove spaces to keep FASTA id clean
    return rid.replace(" ", "_")


def process_record(rec, min_orf_len: int, write_all_frames: bool) -> List[Tuple[str, str]]:
    """Process one SeqRecord; return list of (header, aa_sequence) strings for FASTA.
    """
    dna: Seq = rec.seq
    record_id = safe_record_id(rec)

    frames = generate_six_frames(dna)

    frames_orfs: List[Tuple[str, List[Tuple[int, int, Seq]]]] = []
    for frame_label, aa in frames:
        orfs = find_orfs_in_aa(aa, min_len=min_orf_len)
        frames_orfs.append((frame_label, orfs))

    best = select_best_frame(frames_orfs)

    outputs: List[Tuple[str, str]] = []
    # Write ORFs
    for frame_label, orfs in frames_orfs:
        for start, end, orf_seq in orfs:
            is_best = best is not None and frame_label == best[0] and orf_seq == best[1][2]
            header = fasta_header(record_id, frame_label, start, end, is_best)
            outputs.append((header, str(orf_seq)))

    # Optionally include raw full-frame translations (for traceability)
    if write_all_frames:
        for frame_label, aa in frames:
            header = f">{record_id} frame={frame_label} full_translation"
            outputs.append((header, str(aa)))

    # If no ORFs found, still include the longest full translation per frame to aid debugging
    if not outputs:
        for frame_label, aa in frames:
            header = f">{record_id} frame={frame_label} no_orf_minlen"
            outputs.append((header, str(aa)))

    return outputs


def write_fasta(output_path: Path, entries: Iterable[Tuple[str, str]]) -> None:
    with output_path.open("w", encoding="utf-8") as fh:
        for header, seq in entries:
            fh.write(header + "\n")
            # Wrap to 60 chars per FASTA line
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")


def main() -> int:
    args = parse_args()
    out_path = Path(args.output)

    all_entries: List[Tuple[str, str]] = []
    input_paths = [Path(p) for p in args.input]
    for path in input_paths:
        if not path.exists():
            print(f"[WARN] Input no existe: {path}", file=sys.stderr)
            continue
        try:
            for rec in SeqIO.parse(str(path), "genbank"):
                # Only process nucleotide sequences; guard if input has aa
                if rec.seq and set(str(rec.seq.upper())) <= set("ACGTUNacgtun") | {"N", "n"}:
                    entries = process_record(
                        rec,
                        min_orf_len=args.min_orf_len,
                        write_all_frames=args.write_all_frames,
                    )
                    all_entries.extend(entries)
                else:
                    print(f"[INFO] Omitiendo registro no nucleotídico: {safe_record_id(rec)}", file=sys.stderr)
        except Exception as exc:
            print(f"[ERROR] Al leer {path}: {exc}", file=sys.stderr)

    if not all_entries:
        print("[ERROR] No se generaron secuencias. Verifique inputs/formato.", file=sys.stderr)
        return 1

    write_fasta(out_path, all_entries)
    print(f"[OK] FASTA escrito en: {out_path} (entradas: {len(all_entries)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


