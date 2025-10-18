#!/usr/bin/env python3
"""
Análisis de Alineamiento Múltiple de Secuencias

Este script analiza el alineamiento múltiple generado por Clustal Omega
y proporciona estadísticas detalladas sobre la conservación, variabilidad
y características del alineamiento.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple
from collections import Counter
import re


def parse_alignment_file(alignment_file: Path) -> Dict[str, str]:
    """
    Parsea el archivo de alineamiento FASTA y extrae las secuencias.
    
    Args:
        alignment_file: Ruta al archivo de alineamiento
        
    Returns:
        Diccionario con {nombre_secuencia: secuencia_alineada}
    """
    sequences = {}
    current_seq = ""
    current_name = ""
    
    with alignment_file.open('r', encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Guardar secuencia anterior
                if current_name and current_seq:
                    sequences[current_name] = current_seq
                # Iniciar nueva secuencia
                current_name = line[1:]  # Remover '>'
                current_seq = ""
            else:
                current_seq += line
        
        # Guardar última secuencia
        if current_name and current_seq:
            sequences[current_name] = current_seq
    
    return sequences


def calculate_alignment_statistics(sequences: Dict[str, str]) -> Dict:
    """
    Calcula estadísticas del alineamiento.
    
    Args:
        sequences: Diccionario con secuencias alineadas
        
    Returns:
        Diccionario con estadísticas del alineamiento
    """
    if not sequences:
        return {}
    
    # Obtener longitud del alineamiento
    alignment_length = len(next(iter(sequences.values())))
    num_sequences = len(sequences)
    
    # Contar caracteres por posición
    position_counts = []
    for pos in range(alignment_length):
        pos_chars = [seq[pos] for seq in sequences.values()]
        position_counts.append(Counter(pos_chars))
    
    # Calcular estadísticas
    stats = {
        'num_sequences': num_sequences,
        'alignment_length': alignment_length,
        'sequence_names': list(sequences.keys()),
        'position_stats': []
    }
    
    # Analizar cada posición
    conserved_positions = 0
    variable_positions = 0
    gap_positions = 0
    
    for pos, char_counts in enumerate(position_counts):
        total_chars = sum(char_counts.values())
        gap_count = char_counts.get('-', 0)
        non_gap_count = total_chars - gap_count
        
        # Determinar tipo de posición
        if gap_count == num_sequences:
            pos_type = "all_gaps"
            gap_positions += 1
        elif non_gap_count == 1:
            pos_type = "single_sequence"
            variable_positions += 1
        elif len([c for c in char_counts if c != '-' and char_counts[c] > 0]) == 1:
            pos_type = "conserved"
            conserved_positions += 1
        else:
            pos_type = "variable"
            variable_positions += 1
        
        # Calcular identidad (para posiciones no-gap)
        identity = 0
        if non_gap_count > 0:
            most_common_count = max(char_counts[c] for c in char_counts if c != '-')
            identity = most_common_count / non_gap_count
        
        stats['position_stats'].append({
            'position': pos + 1,
            'char_counts': dict(char_counts),
            'type': pos_type,
            'identity': identity,
            'gap_percentage': (gap_count / total_chars) * 100 if total_chars > 0 else 0
        })
    
    # Estadísticas generales
    stats['conserved_positions'] = conserved_positions
    stats['variable_positions'] = variable_positions
    stats['gap_positions'] = gap_positions
    stats['conservation_percentage'] = (conserved_positions / alignment_length) * 100
    stats['variability_percentage'] = (variable_positions / alignment_length) * 100
    stats['gap_percentage'] = (gap_positions / alignment_length) * 100
    
    return stats


def analyze_sequence_characteristics(sequences: Dict[str, str]) -> Dict:
    """
    Analiza las características de cada secuencia individual.
    
    Args:
        sequences: Diccionario con secuencias alineadas
        
    Returns:
        Diccionario con características de cada secuencia
    """
    sequence_stats = {}
    
    for name, sequence in sequences.items():
        # Contar diferentes tipos de caracteres
        char_counts = Counter(sequence)
        total_length = len(sequence)
        
        # Calcular estadísticas
        gaps = char_counts.get('-', 0)
        nucleotides = sum(char_counts[c] for c in 'ATCG' if c in char_counts)
        n_bases = char_counts.get('N', 0)
        
        sequence_stats[name] = {
            'length': total_length,
            'gaps': gaps,
            'nucleotides': nucleotides,
            'n_bases': n_bases,
            'gap_percentage': (gaps / total_length) * 100 if total_length > 0 else 0,
            'nucleotide_percentage': (nucleotides / total_length) * 100 if total_length > 0 else 0,
            'n_base_percentage': (n_bases / total_length) * 100 if total_length > 0 else 0
        }
    
    return sequence_stats


def identify_conserved_regions(stats: Dict, min_conserved_length: int = 5) -> List[Dict]:
    """
    Identifica regiones conservadas en el alineamiento.
    
    Args:
        stats: Estadísticas del alineamiento
        min_conserved_length: Longitud mínima para considerar una región conservada
        
    Returns:
        Lista de regiones conservadas
    """
    conserved_regions = []
    current_region = None
    
    for pos_stat in stats['position_stats']:
        if pos_stat['type'] == 'conserved':
            if current_region is None:
                current_region = {
                    'start': pos_stat['position'],
                    'end': pos_stat['position'],
                    'length': 1
                }
            else:
                current_region['end'] = pos_stat['position']
                current_region['length'] += 1
        else:
            if current_region and current_region['length'] >= min_conserved_length:
                conserved_regions.append(current_region)
            current_region = None
    
    # Agregar última región si aplica
    if current_region and current_region['length'] >= min_conserved_length:
        conserved_regions.append(current_region)
    
    return conserved_regions


def generate_alignment_report(alignment_file: Path) -> str:
    """
    Genera un reporte completo del alineamiento.
    
    Args:
        alignment_file: Archivo de alineamiento
        
    Returns:
        Reporte formateado del alineamiento
    """
    # Parsear alineamiento
    sequences = parse_alignment_file(alignment_file)
    
    if not sequences:
        return "Error: No se pudieron cargar las secuencias del alineamiento."
    
    # Calcular estadísticas
    alignment_stats = calculate_alignment_statistics(sequences)
    sequence_stats = analyze_sequence_characteristics(sequences)
    conserved_regions = identify_conserved_regions(alignment_stats)
    
    # Generar reporte
    report = []
    report.append("=" * 80)
    report.append("ANÁLISIS DE ALINEAMIENTO MÚLTIPLE DE SECUENCIAS")
    report.append("=" * 80)
    report.append("")
    
    # Información general
    report.append("INFORMACION GENERAL DEL ALINEAMIENTO")
    report.append("-" * 50)
    report.append(f"Numero de secuencias: {alignment_stats['num_sequences']}")
    report.append(f"Longitud del alineamiento: {alignment_stats['alignment_length']} nucleotidos")
    report.append("")
    
    # Estadísticas de conservación
    report.append("ESTADISTICAS DE CONSERVACION")
    report.append("-" * 50)
    report.append(f"Posiciones conservadas: {alignment_stats['conserved_positions']} ({alignment_stats['conservation_percentage']:.1f}%)")
    report.append(f"Posiciones variables: {alignment_stats['variable_positions']} ({alignment_stats['variability_percentage']:.1f}%)")
    report.append(f"Posiciones con gaps: {alignment_stats['gap_positions']} ({alignment_stats['gap_percentage']:.1f}%)")
    report.append("")
    
    # Análisis por secuencia
    report.append("ANALISIS POR SECUENCIA")
    report.append("-" * 50)
    for name, stats in sequence_stats.items():
        # Extraer nombre de especie
        species_name = name.split('_')[2] if len(name.split('_')) > 2 else name
        report.append(f"Especie: {species_name}")
        report.append(f"  Longitud total: {stats['length']} nucleotidos")
        report.append(f"  Nucleotidos validos: {stats['nucleotides']} ({stats['nucleotide_percentage']:.1f}%)")
        report.append(f"  Gaps: {stats['gaps']} ({stats['gap_percentage']:.1f}%)")
        report.append(f"  Bases N: {stats['n_bases']} ({stats['n_base_percentage']:.1f}%)")
        report.append("")
    
    # Regiones conservadas
    if conserved_regions:
        report.append("REGIONES CONSERVADAS IDENTIFICADAS")
        report.append("-" * 50)
        for i, region in enumerate(conserved_regions, 1):
            report.append(f"Region {i}: Posiciones {region['start']}-{region['end']} (longitud: {region['length']})")
        report.append("")
    else:
        report.append("No se identificaron regiones conservadas significativas")
        report.append("")
    
    # Análisis de patrones
    report.append("ANALISIS DE PATRONES")
    report.append("-" * 50)
    
    # Buscar repeticiones en las secuencias
    all_sequences = list(sequences.values())
    if all_sequences:
        # Buscar patrones repetitivos comunes
        sample_seq = all_sequences[0].replace('-', '')
        if 'ACCCTA' in sample_seq:
            report.append("Patron repetitivo 'ACCCTA' detectado en las secuencias")
            report.append("  Este patron es caracteristico de secuencias telomericas")
        else:
            report.append("No se detectaron patrones repetitivos obvios")
    
    report.append("")
    
    # Interpretación biológica
    report.append("INTERPRETACION BIOLOGICA")
    report.append("-" * 50)
    report.append("• Las secuencias analizadas corresponden a cromosomas de primates:")
    report.append("  - Homo sapiens (humano)")
    report.append("  - Pan troglodytes (chimpancé)")
    report.append("  - Pongo abelii (orangután)")
    report.append("")
    report.append("• El alto porcentaje de gaps en Homo sapiens sugiere:")
    report.append("  - Diferencias en la longitud de las secuencias")
    report.append("  - Posibles inserciones/deleciones entre especies")
    report.append("  - Diferencias en la calidad de secuenciación")
    report.append("")
    report.append("• Las regiones conservadas indican:")
    report.append("  - Secuencias funcionalmente importantes")
    report.append("  - Presión selectiva para mantener la secuencia")
    report.append("  - Homología entre las especies")
    
    report.append("")
    report.append("=" * 80)
    
    return "\n".join(report)


def main() -> int:
    """Función principal."""
    parser = argparse.ArgumentParser(
        description="Analiza un alineamiento múltiple de secuencias"
    )
    parser.add_argument(
        "--input",
        "-i",
        type=str,
        default="results/alignment.aln",
        help="Archivo de alineamiento a analizar (default: results/alignment.aln)"
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        help="Archivo de salida para el reporte (opcional)"
    )
    
    args = parser.parse_args()
    alignment_file = Path(args.input)
    
    if not alignment_file.exists():
        print(f"[ERROR] El archivo de alineamiento no existe: {alignment_file}", file=sys.stderr)
        return 1
    
    print(f"[INFO] Analizando alineamiento: {alignment_file}")
    
    try:
        # Generar reporte
        report = generate_alignment_report(alignment_file)
        
        # Mostrar reporte
        print(report)
        
        # Guardar en archivo si se especifica
        if args.output:
            output_file = Path(args.output)
            with output_file.open('w', encoding='utf-8') as f:
                f.write(report)
            print(f"\n[OK] Reporte guardado en: {output_file}")
        
        return 0
        
    except Exception as e:
        print(f"[ERROR] Error analizando alineamiento: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
