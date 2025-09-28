#!/usr/bin/env python3
"""
Ejercicio 2 – Análisis BLAST de secuencias de aminoácidos

Características:
- Lee un archivo FASTA con secuencias de aminoácidos (generado por ex1.py)
- Realiza búsquedas BLASTp contra la base de datos de proteínas de NCBI
- Analiza los resultados y genera un reporte con las mejores coincidencias
- Identifica proteínas conocidas y sus funciones

Uso:
    python ex2.py --input orfs.fasta --output blast_results.txt [--max_hits 10] [--evalue 0.001]

Notas:
- Requiere conexión a internet para acceder a NCBI BLAST
- Utiliza Biopython para realizar búsquedas BLAST online
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import List, Tuple, Dict, Any

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Realiza búsquedas BLASTp de secuencias de aminoácidos "
            "contra la base de datos de proteínas de NCBI."
        )
    )
    parser.add_argument(
        "--input",
        "-i",
        required=True,
        help="Archivo FASTA con secuencias de aminoácidos (generado por ex1.py)",
    )
    parser.add_argument(
        "--output",
        "-o",
        required=True,
        help="Archivo de salida con resultados BLAST",
    )
    parser.add_argument(
        "--max_hits",
        type=int,
        default=10,
        help="Número máximo de hits a reportar por secuencia (default: 10)",
    )
    parser.add_argument(
        "--evalue",
        type=float,
        default=0.001,
        help="Umbral de valor E para filtrar resultados (default: 0.001)",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=1.0,
        help="Delay entre búsquedas BLAST en segundos (default: 1.0)",
    )
    return parser.parse_args()


def format_sequence_for_blast(seq_record) -> str:
    """Formatea una secuencia para BLAST en formato FASTA."""
    return f">{seq_record.id}\n{str(seq_record.seq)}"


def perform_blast_search(sequence: str, program: str = "blastp", database: str = "nr") -> str:
    """
    Realiza una búsqueda BLAST online.
    
    Args:
        sequence: Secuencia en formato FASTA
        program: Programa BLAST a usar (blastp, blastn, etc.)
        database: Base de datos a consultar (nr, nt, etc.)
    
    Returns:
        XML con resultados BLAST
    """
    print(f"Realizando búsqueda BLAST {program} contra {database}...")
    
    try:
        result_handle = NCBIWWW.qblast(program, database, sequence)
        return result_handle.read()
    except Exception as e:
        print(f"Error en búsqueda BLAST: {e}", file=sys.stderr)
        return ""


def parse_blast_results(xml_data: str, max_hits: int = 10, evalue_threshold: float = 0.001) -> List[Dict[str, Any]]:
    """
    Parsea los resultados XML de BLAST y extrae información relevante.
    
    Args:
        xml_data: XML con resultados BLAST
        max_hits: Número máximo de hits a retornar
        evalue_threshold: Umbral de valor E
    
    Returns:
        Lista de diccionarios con información de hits
    """
    if not xml_data:
        return []
    
    hits = []
    
    try:
        from io import StringIO
        blast_record = NCBIXML.read(StringIO(xml_data))
        
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect <= evalue_threshold:
                    hit_info = {
                        'title': alignment.title,
                        'length': alignment.length,
                        'evalue': hsp.expect,
                        'identity': hsp.identities,
                        'positive': hsp.positives,
                        'gaps': hsp.gaps,
                        'query_start': hsp.query_start,
                        'query_end': hsp.query_end,
                        'subject_start': hsp.sbjct_start,
                        'subject_end': hsp.sbjct_end,
                        'query_seq': hsp.query,
                        'match_seq': hsp.match,
                        'subject_seq': hsp.sbjct,
                        'score': hsp.score,
                        'bits': hsp.bits
                    }
                    hits.append(hit_info)
                    
                    if len(hits) >= max_hits:
                        break
            
            if len(hits) >= max_hits:
                break
                
    except Exception as e:
        print(f"Error parseando resultados BLAST: {e}", file=sys.stderr)
    
    return hits


def format_blast_results(sequence_id: str, hits: List[Dict[str, Any]]) -> str:
    """
    Formatea los resultados BLAST para el reporte final.
    
    Args:
        sequence_id: ID de la secuencia consultada
        hits: Lista de hits de BLAST
    
    Returns:
        String formateado con los resultados
    """
    if not hits:
        return f"\n=== BLAST Results for {sequence_id} ===\nNo significant hits found.\n"
    
    output = f"\n=== BLAST Results for {sequence_id} ===\n"
    output += f"Found {len(hits)} significant hits:\n\n"
    
    for i, hit in enumerate(hits, 1):
        # Extraer información básica del título
        title_parts = hit['title'].split('|')
        if len(title_parts) >= 4:
            protein_id = title_parts[3] if len(title_parts) > 3 else "Unknown"
            description = title_parts[-1].strip() if title_parts else "No description"
        else:
            protein_id = "Unknown"
            description = hit['title']
        
        output += f"Hit #{i}:\n"
        output += f"  Protein ID: {protein_id}\n"
        output += f"  Description: {description}\n"
        output += f"  Length: {hit['length']} aa\n"
        output += f"  E-value: {hit['evalue']:.2e}\n"
        output += f"  Score: {hit['score']} (Bits: {hit['bits']:.1f})\n"
        output += f"  Identity: {hit['identity']} / {len(hit['query_seq'])} ({hit['identity']/len(hit['query_seq'])*100:.1f}%)\n"
        output += f"  Positive: {hit['positive']} / {len(hit['query_seq'])} ({hit['positive']/len(hit['query_seq'])*100:.1f}%)\n"
        output += f"  Query range: {hit['query_start']}-{hit['query_end']}\n"
        output += f"  Subject range: {hit['subject_start']}-{hit['subject_end']}\n"
        
        # Mostrar alineamiento (primeros 60 caracteres)
        query_preview = hit['query_seq'][:60]
        match_preview = hit['match_seq'][:60]
        subject_preview = hit['subject_seq'][:60]
        
        output += f"  Alignment preview:\n"
        output += f"    Query:  {query_preview}\n"
        output += f"    Match:  {match_preview}\n"
        output += f"    Subject: {subject_preview}\n"
        output += "\n"
    
    return output


def process_sequences(input_file: Path, max_hits: int, evalue_threshold: float, delay: float) -> str:
    """
    Procesa todas las secuencias del archivo FASTA y realiza búsquedas BLAST.
    
    Args:
        input_file: Archivo FASTA de entrada
        max_hits: Número máximo de hits por secuencia
        evalue_threshold: Umbral de valor E
        delay: Delay entre búsquedas
    
    Returns:
        String con todos los resultados
    """
    all_results = []
    
    try:
        sequences = list(SeqIO.parse(str(input_file), "fasta"))
        print(f"Procesando {len(sequences)} secuencias...")
        
        for i, seq_record in enumerate(sequences, 1):
            print(f"\nProcesando secuencia {i}/{len(sequences)}: {seq_record.id}")
            
            # Formatear secuencia para BLAST
            fasta_sequence = format_sequence_for_blast(seq_record)
            
            # Realizar búsqueda BLAST
            xml_results = perform_blast_search(fasta_sequence)
            
            if xml_results:
                # Parsear resultados
                hits = parse_blast_results(xml_results, max_hits, evalue_threshold)
                
                # Formatear resultados
                formatted_results = format_blast_results(seq_record.id, hits)
                all_results.append(formatted_results)
            else:
                all_results.append(f"\n=== BLAST Results for {seq_record.id} ===\nError: No se pudieron obtener resultados BLAST.\n")
            
            # Delay entre búsquedas para no sobrecargar NCBI
            if i < len(sequences):
                print(f"Esperando {delay} segundos antes de la siguiente búsqueda...")
                time.sleep(delay)
                
    except Exception as e:
        print(f"Error procesando archivo: {e}", file=sys.stderr)
        return f"Error: {e}\n"
    
    return "\n".join(all_results)


def write_results(output_file: Path, results: str) -> None:
    """Escribe los resultados en el archivo de salida."""
    with output_file.open("w", encoding="utf-8") as fh:
        fh.write("BLAST Analysis Results\n")
        fh.write("=" * 50 + "\n")
        fh.write(f"Generated on: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        fh.write("\n")
        fh.write(results)


def main() -> int:
    args = parse_args()
    input_file = Path(args.input)
    output_file = Path(args.output)
    
    if not input_file.exists():
        print(f"[ERROR] Archivo de entrada no existe: {input_file}", file=sys.stderr)
        return 1
    
    print(f"[INFO] Iniciando análisis BLAST...")
    print(f"[INFO] Archivo de entrada: {input_file}")
    print(f"[INFO] Archivo de salida: {output_file}")
    print(f"[INFO] Máximo hits por secuencia: {args.max_hits}")
    print(f"[INFO] Umbral de E-value: {args.evalue}")
    print(f"[INFO] Delay entre búsquedas: {args.delay}s")
    
    # Procesar secuencias
    results = process_sequences(input_file, args.max_hits, args.evalue, args.delay)
    
    # Escribir resultados
    write_results(output_file, results)
    
    print(f"\n[OK] Análisis BLAST completado. Resultados guardados en: {output_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
