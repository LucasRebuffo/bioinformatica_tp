#!/usr/bin/env python3
"""
Script principal para ejecutar el análisis completo de bioinformática.

Este script automatiza la ejecución de ambos ejercicios en secuencia:
1. Ejercicio 1: Procesamiento de secuencias GenBank
2. Ejercicio 2: Análisis BLAST de proteínas

Uso:
    python run_analysis.py --input data/NM_022555.gb --output results/
"""

import argparse
import subprocess
import sys
from pathlib import Path


def run_command(cmd: list, description: str) -> bool:
    """Ejecuta un comando y maneja errores."""
    print(f"\n{'='*60}")
    print(f"Ejecutando: {description}")
    print(f"Comando: {' '.join(cmd)}")
    print(f"{'='*60}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("✓ Comando ejecutado exitosamente")
        if result.stdout:
            print("Salida:", result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        print(f"✗ Error ejecutando comando: {e}")
        if e.stderr:
            print("Error:", e.stderr)
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Ejecuta el análisis completo de bioinformática"
    )
    parser.add_argument(
        "--input",
        "-i",
        required=True,
        help="Archivo GenBank de entrada"
    )
    parser.add_argument(
        "--output-dir",
        "-o",
        default="results",
        help="Directorio de salida (default: results)"
    )
    parser.add_argument(
        "--min-orf-len",
        type=int,
        default=30,
        help="Longitud mínima de ORF (default: 30)"
    )
    parser.add_argument(
        "--max-hits",
        type=int,
        default=10,
        help="Máximo número de hits BLAST (default: 10)"
    )
    parser.add_argument(
        "--evalue",
        type=float,
        default=0.001,
        help="Umbral de E-value para BLAST (default: 0.001)"
    )
    parser.add_argument(
        "--skip-blast",
        action="store_true",
        help="Saltar análisis BLAST (solo ejecutar ejercicio 1)"
    )
    
    args = parser.parse_args()
    
    # Verificar que el archivo de entrada existe
    input_file = Path(args.input)
    if not input_file.exists():
        print(f"Error: El archivo de entrada no existe: {input_file}")
        return 1
    
    # Crear directorio de salida si no existe
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Archivos de salida
    orfs_file = output_dir / "orfs.fasta"
    blast_file = output_dir / "blast_analysis.txt"
    
    print("🧬 Iniciando análisis de bioinformática")
    print(f"📁 Archivo de entrada: {input_file}")
    print(f"📁 Directorio de salida: {output_dir}")
    
    # Ejercicio 1: Procesamiento de secuencias
    cmd1 = [
        "python", "scripts/ex1.py",
        "--input", str(input_file),
        "--output", str(orfs_file),
        "--min_orf_len", str(args.min_orf_len)
    ]
    
    if not run_command(cmd1, "Ejercicio 1: Procesamiento de secuencias GenBank"):
        print("❌ Falló el Ejercicio 1")
        return 1
    
    print(f"✅ Ejercicio 1 completado. ORFs guardados en: {orfs_file}")
    
    # Ejercicio 2: Análisis BLAST (opcional)
    if not args.skip_blast:
        cmd2 = [
            "python", "scripts/ex2.py",
            "--input", str(orfs_file),
            "--output", str(blast_file),
            "--max_hits", str(args.max_hits),
            "--evalue", str(args.evalue)
        ]
        
        if not run_command(cmd2, "Ejercicio 2: Análisis BLAST"):
            print("❌ Falló el Ejercicio 2")
            return 1
        
        print(f"✅ Ejercicio 2 completado. Análisis BLAST guardado en: {blast_file}")
    else:
        print("⏭️  Saltando análisis BLAST (--skip-blast especificado)")
    
    print("\n🎉 Análisis completado exitosamente!")
    print(f"📊 Resultados disponibles en: {output_dir}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
