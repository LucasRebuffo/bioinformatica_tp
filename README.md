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

## Características Técnicas

### Mejoras Implementadas

1. **Manejo de Codones Parciales**: 
   - Los scripts recortan automáticamente las secuencias a múltiplos de 3
   - Elimina advertencias de Biopython sobre codones incompletos

2. **Interfaz en Español**:
   - Todos los comentarios del código están en español
   - Los resultados BLAST se generan en español
   - Mensajes de error y progreso en español

3. **Arquitectura Modular**:
   - Funciones especializadas para cada tarea
   - Código reutilizable y mantenible
   - Separación clara de responsabilidades

4. **Manejo Robusto de Errores**:
   - Validación de archivos de entrada
   - Manejo de errores de conexión BLAST
   - Mensajes informativos detallados

5. **Optimización de Recursos**:
   - Delays configurables entre búsquedas BLAST
   - Filtrado eficiente de resultados
   - Procesamiento por lotes de secuencias

## Notas Importantes

- **Conexión a Internet:** El Ejercicio 2 requiere conexión a internet para acceder a los servidores BLAST de NCBI
- **Rate Limiting:** Se incluye un delay entre búsquedas BLAST para no sobrecargar los servidores de NCBI
- **Formato de Salida:** Los resultados BLAST incluyen información detallada de alineamientos, valores E, identidad y descripciones de proteínas
- **Manejo de Errores:** Ambos scripts incluyen manejo robusto de errores y mensajes informativos

## Dependencias

- `biopython>=1.83,<2`: Para manipulación de secuencias y BLAST
- `python>=3.9`: Versión mínima de Python requerida

## Arquitectura de los Scripts

### Ejercicio 1 (scripts/ex1.py) - Arquitectura Modular

El script está estructurado en funciones especializadas:

1. **`parse_args()`**: Manejo de argumentos de línea de comandos
2. **`generate_six_frames(dna)`**: Genera los 6 marcos de lectura (3 directos + 3 reverso-complementarios)
   - Recorta secuencias a múltiplos de 3 para evitar advertencias de codones parciales
   - Traduce cada marco a aminoácidos
3. **`find_orfs_in_aa(aa_seq, min_len)`**: Detecta ORFs en secuencias de aminoácidos
   - Busca codones de inicio 'M' y parada '*'
   - Filtra por longitud mínima
4. **`select_best_frame(frames_orfs)`**: Selecciona el mejor marco por longitud de ORF
5. **`process_record(rec, min_orf_len, write_all_frames)`**: Procesa un registro GenBank completo
6. **`write_fasta(output_path, entries)`**: Escribe resultados en formato FASTA

**Flujo de datos:**
```
GenBank → 6 Marcos → Traducción → Detección ORFs → Selección Mejor → FASTA
```

### Ejercicio 2 (scripts/ex2.py) - Arquitectura de Análisis BLAST

El script implementa un pipeline de análisis BLAST:

1. **`parse_args()`**: Configuración de parámetros BLAST
2. **`format_sequence_for_blast(seq_record)`**: Formatea secuencias para BLAST
3. **`perform_blast_search(sequence, program, database)`**: Realiza búsquedas BLAST online
   - Utiliza NCBIWWW.qblast() para consultas remotas
   - Maneja errores de conexión
4. **`parse_blast_results(xml_data, max_hits, evalue_threshold)`**: Parsea resultados XML
   - Extrae información de alineamientos
   - Filtra por valor E y número de hits
5. **`format_blast_results(sequence_id, hits)`**: Formatea resultados para reporte
6. **`process_sequences(input_file, max_hits, evalue_threshold, delay)`**: Pipeline principal
   - Procesa múltiples secuencias
   - Implementa delays para no sobrecargar NCBI
7. **`write_results(output_file, results)`**: Genera reporte final

**Flujo de datos:**
```
FASTA → BLAST Online → XML → Parseo → Filtrado → Reporte Español
```

## Estructura de Resultados

### Ejercicio 1 - Salida FASTA
```
>NM_022555.4 frame=+3 aa_start=35 aa_end=301 best_frame=true
MVCLKLPGGSSLAALTVTLMVLSSRLAFAGDTRPRFLELRKSECHFFNGTERVRYLDRYF...
```

### Ejercicio 2 - Análisis BLAST (Resultados en Español)
```
=== Resultados BLAST para NM_022555.4 ===
Se encontraron 5 coincidencias significativas:

Coincidencia #1:
  ID de Proteína: NP_071852.1
  Descripción: HLA class II histocompatibility antigen, DR beta 3 chain
  Longitud: 266 aa
  Valor E: 0.0
  Puntuación: 266 (Bits: 266.0)
  Identidad: 266 / 266 (100.0%)
  ...
```

## Contribuciones

Este proyecto forma parte de un trabajo práctico de bioinformática. Para modificaciones o mejoras, contactar al autor del proyecto.


