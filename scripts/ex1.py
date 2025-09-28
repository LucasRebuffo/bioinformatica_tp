#!/usr/bin/env python3
"""
Ejercicio 1 – Procesamiento de secuencias de GenBank (registros de ARNm NM_...)

Características:
- Lee uno o más archivos GenBank que contienen secuencias de referencia de ARNm (NM_...)
- Para cada registro, genera los 6 marcos de lectura (3 directos, 3 reverso-complementarios)
- Traduce los marcos a secuencias de aminoácidos
- Detecta marcos de lectura abiertos (ORFs) usando el codón de inicio 'M' y de parada '*'
- Selecciona el mejor marco como aquel que contiene el ORF más largo (configurable)
- Escribe las secuencias de aminoácidos en un archivo de salida FASTA, incluyendo la anotación del mejor marco

Uso:
    python ex1.py --input NMxxxx.gbk [NMyyyy.gbk ...] --output output.fasta

Si un archivo GenBank contiene múltiples registros, cada uno será procesado.

Notas:
- Requiere Biopython.
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
    """Devuelve los seis marcos de lectura como (etiqueta_marco, sec_aa) sin eliminar codones de parada.

    Los marcos se etiquetan: +1, +2, +3, -1, -2, -3
    """
    frames: List[Tuple[str, Seq]] = []
    # Marcos de lectura directos
    for offset in range(3):
        frame_nt = dna[offset:]
        # Recortar a múltiplo de 3 para evitar advertencia de codón parcial
        frame_nt = frame_nt[:len(frame_nt) - (len(frame_nt) % 3)]
        aa = frame_nt.translate(to_stop=False)
        frames.append((f"+{offset+1}", aa))
    # Marcos de lectura reverso-complementarios
    rc = dna.reverse_complement()
    for offset in range(3):
        frame_nt = rc[offset:]
        # Recortar a múltiplo de 3 para evitar advertencia de codón parcial
        frame_nt = frame_nt[:len(frame_nt) - (len(frame_nt) % 3)]
        aa = frame_nt.translate(to_stop=False)
        frames.append((f"-{offset+1}", aa))
    return frames


def find_orfs_in_aa(aa_seq: Seq, min_len: int) -> List[Tuple[int, int, Seq]]:
    """Encuentra ORFs en una secuencia de aminoácidos.

    Devuelve una lista de (indice_inicio, indice_fin, sec_orf), donde los índices son base 0 y el final es exclusivo.
    Definición de ORF: comienza en 'M', termina en el siguiente '*' o al final de la secuencia. Debe ser >= min_len aa.
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
    # Si la secuencia termina en un ORF sin codón de parada terminal, considerarlo también
    if start_pos is not None:
        if len(aa_seq) - start_pos >= min_len:
            orf_seq = aa_seq[start_pos:]
            orfs.append((start_pos, len(aa_seq), orf_seq))
    return orfs


def select_best_frame(frames_orfs: List[Tuple[str, List[Tuple[int, int, Seq]]]]) -> Tuple[str, Tuple[int, int, Seq]] | None:
    """Selecciona el mejor marco por la longitud del ORF más largo.

    Devuelve (etiqueta_marco, (inicio, fin, sec_orf)) o None si no se encuentran ORFs.
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
    # Preferir accession o id
    rid = getattr(rec, "id", None) or getattr(rec, "name", None) or "unknown"
    # Eliminar espacios para mantener limpio el id de FASTA
    return rid.replace(" ", "_")


def process_record(rec, min_orf_len: int, write_all_frames: bool) -> List[Tuple[str, str]]:
    """Procesa un SeqRecord; devuelve una lista de tuplas (cabecera, secuencia_aa) para el FASTA.
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
    # Escribir ORFs
    for frame_label, orfs in frames_orfs:
        for start, end, orf_seq in orfs:
            is_best = best is not None and frame_label == best[0] and orf_seq == best[1][2]
            header = fasta_header(record_id, frame_label, start, end, is_best)
            outputs.append((header, str(orf_seq)))

    # Opcionalmente, incluir traducciones completas de los marcos (para trazabilidad)
    if write_all_frames:
        for frame_label, aa in frames:
            header = f">{record_id} frame={frame_label} full_translation"
            outputs.append((header, str(aa)))

    # Si no se encuentran ORFs, incluir la traducción completa más larga por marco para ayudar a depurar
    if not outputs:
        for frame_label, aa in frames:
            header = f">{record_id} frame={frame_label} no_orf_minlen"
            outputs.append((header, str(aa)))

    return outputs


def write_fasta(output_path: Path, entries: Iterable[Tuple[str, str]]) -> None:
    with output_path.open("w", encoding="utf-8") as fh:
        for header, seq in entries:
            fh.write(header + "\n")
            # Ajustar a 60 caracteres por línea de FASTA
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
                # Solo procesar secuencias de nucleótidos; protección si el input tiene aa
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


