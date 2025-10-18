#!/usr/bin/env python3
"""
Ejercicio 3 – Alineamiento múltiple de secuencias con Clustal Omega

Características:
- Lee archivos FASTA de secuencias de nucleótidos desde un directorio
- Combina las secuencias en un archivo FASTA limpio
- Realiza alineamiento múltiple usando Clustal Omega
- Genera archivo de alineamiento en formato FASTA

Uso:
    python ex3.py --input ../data/secuencias --output secuencias_alineadas.aln [--combined secuencias_combinadas.fa]

Notas:
- Requiere Biopython y Clustal Omega instalado
- Las secuencias se limpian automáticamente (solo nucleótidos válidos)
- Los IDs de secuencias se normalizan (espacios reemplazados por guiones bajos)
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
from pathlib import Path
from typing import List

from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Seq import Seq


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Realiza alineamiento múltiple de secuencias FASTA "
            "usando Clustal Omega."
        )
    )
    parser.add_argument(
        "--input",
        "-i",
        type=str,
        default="./data/secuencias",
        help="Directorio que contiene archivos FASTA (default: ../data/secuencias)",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="./results/secuencias_alineadas.aln",
        help="Archivo de salida con alineamiento (default: secuencias_alineadas.aln)",
    )
    parser.add_argument(
        "--combined",
        type=str,
        default="./results/secuencias_combinadas.fa",
        help="Archivo FASTA combinado temporal (default: secuencias_combinadas.fa)",
    )
    parser.add_argument(
        "--clustalo",
        type=str,
        help="Ruta al ejecutable de Clustal Omega (auto-detecta si no se especifica)",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Habilitar salida detallada de Clustal Omega",
    )
    return parser.parse_args()


def clean_sequence_id(record_id: str) -> str:
    """Limpia y normaliza el ID de una secuencia."""
    return record_id.strip().replace(" ", "_")


def clean_sequence_nucleotides(sequence: str) -> str:
    """Limpia una secuencia de nucleótidos, manteniendo solo caracteres válidos."""
    return ''.join(ch for ch in str(sequence).upper() if ch in "ACGTU")


def find_fasta_files(input_dir: str) -> List[str]:
    """Encuentra archivos FASTA en el directorio especificado."""
    pattern = os.path.join(input_dir, "*.fasta")
    files = glob.glob(pattern)
    
    if not files:
        raise FileNotFoundError(f"No se encontraron archivos FASTA en {input_dir}")
    
    return files


def combine_fasta_files(input_files: List[str], output_file: str) -> int:
    """
    Combina archivos FASTA en un solo archivo limpio.
    
    Args:
        input_files: Lista de archivos FASTA de entrada
        output_file: Archivo FASTA de salida combinado
    
    Returns:
        Número de secuencias procesadas
    """
    sequences_count = 0
    
    with open(output_file, "w", encoding="ascii", newline="\n") as output_handle:
        for file_path in input_files:
            try:
                for record in SeqIO.parse(file_path, "fasta"):
                    # Limpiar ID y secuencia
                    clean_id = clean_sequence_id(record.id)
                    clean_seq = clean_sequence_nucleotides(record.seq)
                    
                    # Crear nuevo registro limpio
                    record.id = clean_id
                    record.description = ""
                    record.seq = Seq(clean_seq)
                    
                    SeqIO.write(record, output_handle, "fasta")
                    sequences_count += 1
                    
            except Exception as e:
                print(f"[WARN] Error procesando {file_path}: {e}", file=sys.stderr)
                continue
    
    return sequences_count


def find_clustal_omega_executable(clustalo_path: str | None) -> str:
    """Encuentra el ejecutable de Clustal Omega."""
    if clustalo_path and os.path.exists(clustalo_path):
        return clustalo_path
    
    # Auto-detectar en el directorio del script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    auto_path = os.path.join(script_dir, "clustal-omega-1.2.2", "clustalo.exe")
    
    if os.path.exists(auto_path):
        return auto_path
    
    # Buscar en PATH
    import shutil
    clustalo_cmd = shutil.which("clustalo")
    if clustalo_cmd:
        return clustalo_cmd
    
    raise FileNotFoundError(
        "No se encontró Clustal Omega. Especifique la ruta con --clustalo o "
        "instale Clustal Omega en el PATH del sistema."
    )


def run_clustal_omega(
    input_file: str, 
    output_file: str, 
    clustalo_exe: str, 
    verbose: bool = False
) -> None:
    """
    Ejecuta Clustal Omega para realizar el alineamiento múltiple.
    
    Args:
        input_file: Archivo FASTA de entrada
        output_file: Archivo de alineamiento de salida
        clustalo_exe: Ruta al ejecutable de Clustal Omega
        verbose: Habilitar salida detallada
    """
    # Configurar el comando
    clustalo_cline = ClustalOmegaCommandline(
        cmd=clustalo_exe,
        infile=input_file,
        outfile=output_file,
        verbose=verbose,
        auto=True,
        force=True
    )
    
    print(f"[INFO] Ejecutando Clustal Omega...")
    print(f"[INFO] Comando: {clustalo_cline}")
    
    try:
        # Ejecutar Clustal Omega
        stdout, stderr = clustalo_cline()
        
        if stdout:
            print(f"[INFO] Salida de Clustal Omega:\n{stdout}")
        if stderr:
            print(f"[INFO] Log de Clustal Omega:\n{stderr}")
            
        print(f"[OK] Clustal Omega ejecutado correctamente")
        
    except Exception as e:
        raise RuntimeError(f"Error ejecutando Clustal Omega: {e}")


def main() -> int:
    """Función principal que orquesta el proceso de alineamiento."""
    args = parse_args()
    
    # Convertir rutas a Path objects
    input_dir = Path(args.input)
    output_file = Path(args.output)
    combined_file = Path(args.combined)
    
    print(f"[INFO] Iniciando alineamiento múltiple con Clustal Omega")
    print(f"[INFO] Directorio de entrada: {input_dir}")
    print(f"[INFO] Archivo de salida: {output_file}")
    print(f"[INFO] Archivo combinado: {combined_file}")
    
    try:
        # Verificar que el directorio de entrada existe
        if not input_dir.exists():
            print(f"[ERROR] El directorio de entrada no existe: {input_dir}", file=sys.stderr)
            return 1
        
        # Encontrar archivos FASTA
        print(f"[INFO] Buscando archivos FASTA en {input_dir}...")
        input_files = find_fasta_files(str(input_dir))
        print(f"[INFO] Encontrados {len(input_files)} archivos FASTA")
        
        # Limpiar archivos previos si existen
        for f in [combined_file, output_file]:
            if f.exists():
                print(f"[INFO] Eliminando archivo previo: {f}")
                f.unlink()
        
        # Combinar archivos FASTA
        print(f"[INFO] Combinando secuencias en {combined_file}...")
        sequences_count = combine_fasta_files(input_files, str(combined_file))
        print(f"[OK] Archivo combinado creado con {sequences_count} secuencias")
        
        if sequences_count == 0:
            print(f"[ERROR] No se procesaron secuencias válidas", file=sys.stderr)
            return 1
        
        # Encontrar ejecutable de Clustal Omega
        print(f"[INFO] Localizando Clustal Omega...")
        clustalo_exe = find_clustal_omega_executable(args.clustalo)
        print(f"[INFO] Usando Clustal Omega: {clustalo_exe}")
        
        # Ejecutar alineamiento
        run_clustal_omega(
            str(combined_file), 
            str(output_file), 
            clustalo_exe, 
            args.verbose
        )
        
        # Verificar que se generó el archivo de alineamiento
        if not output_file.exists():
            print(f"[ERROR] No se generó el archivo de alineamiento: {output_file}", file=sys.stderr)
            return 1
        
        print(f"[OK] Alineamiento completado exitosamente")
        print(f"[OK] Archivo de alineamiento: {output_file}")
        print(f"[INFO] Secuencias procesadas: {sequences_count}")
        
        return 0
        
    except FileNotFoundError as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"[ERROR] Error inesperado: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())