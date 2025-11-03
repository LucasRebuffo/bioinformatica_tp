#!/usr/bin/env python3
"""
Ejercicio 4 – Análisis de reporte BLAST con filtrado por patrón

Características:
- Parsea un reporte de salida de BLAST (formato texto)
- Identifica hits que contengan un patrón específico en su descripción
- Extrae el ACCESSION de los hits identificados
- Obtiene las secuencias completas de GenBank usando Bio.Entrez
- Escribe las secuencias en formato FASTA

Uso:
    python ex4.py --input blast_results.txt --pattern "Homo sapiens" --output filtered_hits.txt [--fasta output.fasta]

Notas:
- Requiere conexión a internet para acceder a GenBank
- Utiliza Biopython para obtener secuencias de NCBI
"""

from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path
from typing import List, Dict, Optional, Tuple

from Bio import Entrez, SeqIO
from io import StringIO


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Parsea un reporte BLAST y filtra hits por patrón en la descripción. "
            "Opcionalmente obtiene las secuencias completas de GenBank."
        )
    )
    parser.add_argument(
        "--input",
        "-i",
        required=True,
        help="Archivo de reporte BLAST (formato txt)",
    )
    parser.add_argument(
        "--pattern",
        "-p",
        required=True,
        help="Patrón a buscar en las descripciones (ej: 'Homo sapiens')",
    )
    parser.add_argument(
        "--output",
        "-o",
        required=True,
        help="Archivo de salida con lista de hits filtrados",
    )
    parser.add_argument(
        "--fasta",
        "-f",
        default=None,
        help="Archivo FASTA opcional para guardar secuencias completas de GenBank",
    )
    parser.add_argument(
        "--email",
        "-e",
        default="lrebuffo@frba.utn.edu.ar",
        help="Email para identificación con NCBI Entrez (requerido por NCBI)",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=0.5,
        help="Delay entre consultas a GenBank en segundos (default: 0.5)",
    )
    return parser.parse_args()


def parse_blast_report(file_path: Path) -> List[Dict[str, str]]:
    """
    Parsea un archivo de reporte BLAST en formato texto.
    
    Args:
        file_path: Ruta al archivo de reporte BLAST
        
    Returns:
        Lista de diccionarios con información de cada hit
    """
    hits = []
    current_hit = {}
    current_query = None
    in_hit = False
    
    try:
        with file_path.open("r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                
                # Detectar inicio de resultados para una secuencia query
                if line.startswith("=== Resultados BLAST para"):
                    match = re.search(r"===\s+Resultados BLAST para\s+(\S+)", line)
                    if match:
                        current_query = match.group(1)
                        in_hit = False
                
                # Detectar inicio de un hit
                elif line.startswith("Coincidencia #"):
                    if current_hit:
                        hits.append(current_hit)
                    # Extraer número del hit (puede ser "Coincidencia #1:" o "Coincidencia #1")
                    hit_num_match = re.search(r'#(\d+)', line)
                    hit_number = hit_num_match.group(1) if hit_num_match else "N/A"
                    current_hit = {
                        "query_id": current_query,
                        "hit_number": hit_number
                    }
                    in_hit = True
                
                # Extraer ID de Proteína
                elif in_hit and line.startswith("ID de Proteína:"):
                    protein_id = line.split(":", 1)[1].strip()
                    current_hit["protein_id"] = protein_id
                    # Intentar extraer accession de la descripción si es Unknown
                    if protein_id == "Unknown":
                        current_hit["accession"] = None
                    else:
                        current_hit["accession"] = protein_id
                
                # Extraer Descripción
                elif in_hit and line.startswith("Descripción:"):
                    description = line.split(":", 1)[1].strip()
                    current_hit["description"] = description
                    
                    # Intentar extraer accession de la descripción si está en formato emb|CAA45118.1| o similar
                    if not current_hit.get("accession"):
                        # Buscar patrones como: emb|CAA45118.1|, ref|NP_001393645.1|, gb|ADO22212.1|, etc.
                        accession_match = re.search(r'(?:emb|ref|gb|dbj)\|([A-Z0-9_\.]+)\|', description)
                        if accession_match:
                            current_hit["accession"] = accession_match.group(1)
                
                # Extraer otras propiedades del hit
                elif in_hit and line.startswith("Longitud:"):
                    length_match = re.search(r'(\d+)\s+aa', line)
                    if length_match:
                        current_hit["length"] = length_match.group(1)
                
                elif in_hit and line.startswith("Valor E:"):
                    evalue_match = re.search(r'([\d\.]+[eE][\+\-]\d+|\d+\.\d+)', line)
                    if evalue_match:
                        current_hit["evalue"] = evalue_match.group(1)
                
                elif in_hit and line.startswith("Puntuación:"):
                    score_match = re.search(r'(\d+\.?\d*)', line)
                    if score_match:
                        current_hit["score"] = score_match.group(1)
            
            # Agregar el último hit
            if current_hit:
                hits.append(current_hit)
                
    except Exception as e:
        print(f"[ERROR] Error parseando archivo BLAST: {e}", file=sys.stderr)
        return []
    
    return hits


def filter_hits_by_pattern(hits: List[Dict[str, str]], pattern: str, case_sensitive: bool = False) -> List[Dict[str, str]]:
    """
    Filtra hits que contengan el patrón en su descripción.
    
    Args:
        hits: Lista de hits
        pattern: Patrón a buscar
        case_sensitive: Si es True, la búsqueda es sensible a mayúsculas/minúsculas
        
    Returns:
        Lista de hits filtrados
    """
    filtered = []
    search_pattern = pattern if case_sensitive else pattern.lower()
    
    for hit in hits:
        description = hit.get("description", "")
        if not case_sensitive:
            description = description.lower()
        
        if search_pattern in description:
            filtered.append(hit)
    
    return filtered


def extract_accessions(hits: List[Dict[str, str]]) -> List[str]:
    """
    Extrae los accessions únicos de la lista de hits.
    
    Args:
        hits: Lista de hits filtrados
        
    Returns:
        Lista de accessions únicos (sin None)
    """
    accessions = []
    for hit in hits:
        acc = hit.get("accession")
        if acc and acc not in accessions:
            accessions.append(acc)
    return accessions


def fetch_sequences_from_genbank(accessions: List[str], email: str, delay: float = 0.5) -> List[Tuple[str, str]]:
    """
    Obtiene secuencias completas de GenBank usando los accessions.
    
    Args:
        accessions: Lista de accessions a buscar
        email: Email para identificación con NCBI
        delay: Delay entre consultas
        
    Returns:
        Lista de tuplas (accession, secuencia_fasta)
    """
    Entrez.email = email
    sequences = []
    
    if not accessions:
        return sequences
    
    print(f"[INFO] Obteniendo {len(accessions)} secuencia(s) de GenBank...")
    
    for i, accession in enumerate(accessions, 1):
        try:
            print(f"[INFO] Obteniendo secuencia {i}/{len(accessions)}: {accession}")
            
            # Buscar en GenBank
            handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
            fasta_data = handle.read()
            handle.close()
            
            if fasta_data:
                sequences.append((accession, fasta_data))
                print(f"[OK] Secuencia obtenida: {accession}")
            else:
                print(f"[WARNING] No se pudo obtener secuencia para {accession}", file=sys.stderr)
            
            # Delay entre consultas para no sobrecargar NCBI
            if i < len(accessions):
                time.sleep(delay)
                
        except Exception as e:
            print(f"[ERROR] Error obteniendo secuencia {accession}: {e}", file=sys.stderr)
            continue
    
    return sequences


def format_filtered_hits(hits: List[Dict[str, str]], pattern: str) -> str:
    """
    Formatea los hits filtrados para el reporte de salida.
    
    Args:
        hits: Lista de hits filtrados
        pattern: Patrón usado para filtrar
        
    Returns:
        String formateado con los resultados
    """
    if not hits:
        return f"No se encontraron hits que coincidan con el patrón '{pattern}'.\n"
    
    output = f"Análisis de Reporte BLAST - Filtrado por Patrón\n"
    output += "=" * 60 + "\n"
    output += f"Patrón buscado: '{pattern}'\n"
    output += f"Total de hits encontrados: {len(hits)}\n\n"
    
    # Agrupar por query
    hits_by_query = {}
    for hit in hits:
        query_id = hit.get("query_id", "Unknown")
        if query_id not in hits_by_query:
            hits_by_query[query_id] = []
        hits_by_query[query_id].append(hit)
    
    for query_id, query_hits in hits_by_query.items():
        output += f"\n=== Query: {query_id} ===\n"
        output += f"Hits encontrados: {len(query_hits)}\n\n"
        
        for hit in query_hits:
            output += f"Hit #{hit.get('hit_number', 'N/A')}:\n"
            output += f"  ID de Proteína: {hit.get('protein_id', 'N/A')}\n"
            output += f"  ACCESSION: {hit.get('accession', 'N/A')}\n"
            output += f"  Descripción: {hit.get('description', 'N/A')}\n"
            output += f"  Longitud: {hit.get('length', 'N/A')} aa\n"
            output += f"  Valor E: {hit.get('evalue', 'N/A')}\n"
            output += f"  Puntuación: {hit.get('score', 'N/A')}\n"
            output += "\n"
    
    return output


def write_fasta_file(file_path: Path, sequences: List[Tuple[str, str]]) -> None:
    """
    Escribe las secuencias en formato FASTA a un archivo.
    
    Args:
        file_path: Ruta al archivo de salida
        sequences: Lista de tuplas (accession, secuencia_fasta)
    """
    if not sequences:
        print("[WARNING] No hay secuencias para escribir al archivo FASTA", file=sys.stderr)
        return
    
    try:
        with file_path.open("w", encoding="utf-8") as f:
            for accession, fasta_data in sequences:
                f.write(fasta_data)
                if not fasta_data.endswith("\n"):
                    f.write("\n")
        
        print(f"[OK] {len(sequences)} secuencia(s) escrita(s) en: {file_path}")
    except Exception as e:
        print(f"[ERROR] Error escribiendo archivo FASTA: {e}", file=sys.stderr)


def main() -> int:
    args = parse_args()
    input_file = Path(args.input)
    output_file = Path(args.output)
    
    if not input_file.exists():
        print(f"[ERROR] Archivo de entrada no existe: {input_file}", file=sys.stderr)
        return 1
    
    print(f"[INFO] Parseando archivo BLAST: {input_file}")
    
    # Parsear archivo BLAST
    all_hits = parse_blast_report(input_file)
    print(f"[INFO] Total de hits parseados: {len(all_hits)}")
    
    # Filtrar hits por patrón
    print(f"[INFO] Filtrando hits con patrón: '{args.pattern}'")
    filtered_hits = filter_hits_by_pattern(all_hits, args.pattern, case_sensitive=False)
    print(f"[INFO] Hits que coinciden con el patrón: {len(filtered_hits)}")
    
    # Formatear y escribir resultados
    results_text = format_filtered_hits(filtered_hits, args.pattern)
    try:
        with output_file.open("w", encoding="utf-8") as f:
            f.write(results_text)
        print(f"[OK] Resultados guardados en: {output_file}")
    except Exception as e:
        print(f"[ERROR] Error escribiendo archivo de salida: {e}", file=sys.stderr)
        return 1
    
    # Si se especificó archivo FASTA, obtener secuencias completas
    if args.fasta:
        fasta_file = Path(args.fasta)
        accessions = extract_accessions(filtered_hits)
        
        if accessions:
            print(f"[INFO] Extrayendo {len(accessions)} accession(s) único(s)...")
            sequences = fetch_sequences_from_genbank(accessions, args.email, args.delay)
            write_fasta_file(fasta_file, sequences)
        else:
            print("[WARNING] No se encontraron accessions para obtener secuencias", file=sys.stderr)
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

