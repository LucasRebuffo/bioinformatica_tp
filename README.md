# Proyecto de Bioinformática - Análisis de Secuencias y BLAST

Este proyecto contiene herramientas de bioinformática para el análisis de secuencias de ARNm, detección de ORFs (Open Reading Frames) y análisis BLAST de proteínas.

## Estructura del Proyecto

```
bioinformatica_tp/
├── scripts/           # Scripts de Python
│   ├── ex1.py        # Ejercicio 1: Procesamiento de secuencias GenBank
│   └── ex2.py        # Ejercicio 2: Análisis BLAST de proteínas
├── data/             # Archivos de datos de entrada
│   ├── NM_022555.gb  # Archivo GenBank de ejemplo
│   └── INS_orfs.fasta # Secuencias de aminoácidos generadas
├── results/          # Archivos de resultados
│   ├── blast_analysis.txt # Resultados de análisis BLAST
│   └── INS_orfs_table.xlsx # Tabla de resultados
├── docs/             # Documentación
│   └── TP Bioinformatica-2025 - Parte 1.pdf
├── run_analysis.py   # Script principal para análisis automatizado
├── requirements.txt   # Dependencias de Python
└── README.md         # Este archivo
```

## Requisitos

- Python 3.9+
- Dependencias: ver `requirements.txt`

### Instalación

```bash
pip install -r requirements.txt
```

## Ejercicio 1 – Procesamiento de Secuencias GenBank

**Script:** `scripts/ex1.py`

Este script lee uno o más archivos GenBank con mRNA(s) de referencia (NM_...), genera las 6 traducciones posibles (3 marcos forward y 3 reverse complement), detecta ORFs y escribe un FASTA con las secuencias de aminoácidos.

### Características
- Genera los 6 marcos de lectura (3 directos, 3 reverso-complementarios)
- Traduce los marcos a secuencias de aminoácidos
- Detecta marcos de lectura abiertos (ORFs) usando el codón de inicio 'M' y de parada '*'
- Selecciona el mejor marco como aquel que contiene el ORF más largo
- Marca como "best_frame=true" al ORF más largo encontrado

### Uso
```bash
python scripts/ex1.py --input data/NMxxxx.gbk --output results/orfs.fasta [--min_orf_len 30] [--write_all_frames]
```

### Parámetros
- `--input/-i`: uno o más archivos GenBank (`.gb`/`.gbk`) con mRNA(s) NM_
- `--output/-o`: archivo FASTA de salida con aminoácidos
- `--min_orf_len`: longitud mínima de ORF en aa (default: 30)
- `--write_all_frames`: incluye traducciones completas de los 6 marcos

### Ejemplo
```bash
python scripts/ex1.py -i data/NM_022555.gb -o results/INS_orfs.fasta --min_orf_len 30
```

## Ejercicio 2 – Análisis BLAST de Proteínas

**Script:** `scripts/ex2.py`

Este script realiza búsquedas BLASTp de secuencias de aminoácidos contra la base de datos de proteínas de NCBI para identificar proteínas similares y sus funciones.

### Características
- Lee archivos FASTA con secuencias de aminoácidos
- Realiza búsquedas BLASTp online contra la base de datos nr de NCBI
- Analiza y filtra resultados por valor E
- Genera reportes detallados con información de alineamientos
- Identifica proteínas conocidas y sus funciones

### Uso
```bash
python scripts/ex2.py --input results/orfs.fasta --output results/blast_results.txt [--max_hits 10] [--evalue 0.001]
```

### Parámetros
- `--input/-i`: archivo FASTA con secuencias de aminoácidos
- `--output/-o`: archivo de salida con resultados BLAST
- `--max_hits`: número máximo de hits a reportar por secuencia (default: 10)
- `--evalue`: umbral de valor E para filtrar resultados (default: 0.001)
- `--delay`: delay entre búsquedas BLAST en segundos (default: 1.0)

### Ejemplo
```bash
python scripts/ex2.py -i results/INS_orfs.fasta -o results/blast_analysis.txt --max_hits 5 --evalue 0.01
```

## Flujo de Trabajo Completo

### Opción 1: Script Automatizado (Recomendado)

Usar el script principal `run_analysis.py` que ejecuta ambos ejercicios automáticamente:

```bash
# Análisis completo
python run_analysis.py -i data/NM_022555.gb -o results/

# Con parámetros personalizados
python run_analysis.py -i data/NM_022555.gb -o results/ --min-orf-len 30 --max-hits 5 --evalue 0.01

# Solo ejercicio 1 (sin BLAST)
python run_analysis.py -i data/NM_022555.gb -o results/ --skip-blast
```

### Opción 2: Ejecución Manual

1. **Preparar datos de entrada:**
   ```bash
   # Descargar archivos GenBank desde NCBI
   # Colocar en carpeta data/
   ```

2. **Ejecutar análisis de secuencias:**
   ```bash
   python scripts/ex1.py -i data/NM_022555.gb -o results/orfs.fasta --min_orf_len 30
   ```

3. **Realizar análisis BLAST:**
   ```bash
   python scripts/ex2.py -i results/orfs.fasta -o results/blast_analysis.txt --max_hits 10
   ```

## Notas Importantes

- **Conexión a Internet:** El Ejercicio 2 requiere conexión a internet para acceder a los servidores BLAST de NCBI
- **Rate Limiting:** Se incluye un delay entre búsquedas BLAST para no sobrecargar los servidores de NCBI
- **Formato de Salida:** Los resultados BLAST incluyen información detallada de alineamientos, valores E, identidad y descripciones de proteínas
- **Manejo de Errores:** Ambos scripts incluyen manejo robusto de errores y mensajes informativos

## Dependencias

- `biopython>=1.83,<2`: Para manipulación de secuencias y BLAST
- `python>=3.9`: Versión mínima de Python requerida

## Estructura de Resultados

### Ejercicio 1 - Salida FASTA
```
>NM_022555.4 frame=+3 aa_start=35 aa_end=301 best_frame=true
MVCLKLPGGSSLAALTVTLMVLSSRLAFAGDTRPRFLELRKSECHFFNGTERVRYLDRYF...
```

### Ejercicio 2 - Análisis BLAST
```
=== BLAST Results for NM_022555.4 ===
Found 5 significant hits:

Hit #1:
  Protein ID: NP_071852.1
  Description: HLA class II histocompatibility antigen, DR beta 3 chain
  Length: 266 aa
  E-value: 0.0
  Score: 266 (Bits: 266.0)
  Identity: 266 / 266 (100.0%)
  ...
```

## Contribuciones

Este proyecto forma parte de un trabajo práctico de bioinformática. Para modificaciones o mejoras, contactar al autor del proyecto.


