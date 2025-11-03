# Proyecto de Bioinformática - Análisis de Secuencias y BLAST

Este proyecto contiene herramientas de bioinformática para el análisis de secuencias de ARNm, detección de ORFs (Open Reading Frames) y análisis BLAST de proteínas.

## Estructura del Proyecto

```
bioinformatica_tp/
├── scripts/                   # Scripts de Python
│   ├── clustal-omega-1.2.2/   # Proyecto Clustal Omega para alinear secuencias
│   ├── ex1.py                 # Ejercicio 1: Procesamiento de secuencias GenBank
│   ├── ex2.py                 # Ejercicio 2: Análisis BLAST de proteínas
│   ├── ex3.py                 # Ejercicio 3: Alineación múltiple de secuencias
│   └── ex4.py                 # Ejercicio 4: Análisis de reporte BLAST con filtrado por patrón
├── data/             # Archivos de datos de entrada
│   ├── secuencias/   # Secuencias fasta para el ejercicio 3
│   ├── NM_022555.gb  # Archivo GenBank de ejemplo
│   └── INS_orfs.fasta # Secuencias de aminoácidos generadas
├── results/          # Archivos de resultados
│   ├── blast_analysis.txt # Resultados de análisis BLAST
│   └── INS_orfs_table.xlsx # Tabla de resultados
├── docs/             # Documentación
│   └── TP Bioinformatica-2025 - Parte 1.pdf
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

## Ejercicio 3 – Alineamiento Múltiple de Secuencias

**Script:** `scripts/ex3.py`

Este script realiza alineamiento múltiple de secuencias de nucleótidos usando Clustal Omega. Combina archivos FASTA de entrada y genera un alineamiento múltiple de alta calidad.

### Características
- Lee archivos FASTA desde un directorio especificado
- Combina múltiples archivos FASTA en un archivo limpio
- Realiza alineamiento múltiple usando Clustal Omega
- Limpia automáticamente secuencias (solo nucleótidos válidos)
- Normaliza IDs de secuencias
- Auto-detecta el ejecutable de Clustal Omega
- Manejo robusto de errores y mensajes informativos

### Uso
```bash
python scripts/ex3.py --input ../data/secuencias --output secuencias_alineadas.aln [--combined secuencias_combinadas.fa] [--clustalo /path/to/clustalo] [--verbose]
```

### Parámetros
- `--input/-i`: directorio que contiene archivos FASTA (default: `../data/secuencias`)
- `--output/-o`: archivo de salida con alineamiento (default: `secuencias_alineadas.aln`)
- `--combined`: archivo FASTA combinado temporal (default: `secuencias_combinadas.fa`)
- `--clustalo`: ruta al ejecutable de Clustal Omega (auto-detecta si no se especifica)
- `--verbose`: habilitar salida detallada de Clustal Omega

### Ejemplo
```bash
python scripts/ex3.py -i ./data/secuencias -o results/alignment.aln --verbose
```

## Ejercicio 4 – Análisis de Reporte BLAST con Filtrado por Patrón

**Script:** `scripts/ex4.py`

Este script parsea un reporte de salida de BLAST (formato texto) e identifica los hits que contienen un patrón específico en su descripción. Además, implementa el punto extra: extrae los ACCESSION de los hits identificados y obtiene las secuencias completas de GenBank en formato FASTA.

### Características
- Parsea archivos de reporte BLAST en formato texto (generados por `ex2.py`)
- Filtra hits cuya descripción contiene un patrón específico (búsqueda case-insensitive)
- Extrae los números de acceso (ACCESSION) de los hits identificados
- Obtiene las secuencias completas de GenBank usando Bio.Entrez
- Escribe las secuencias en formato FASTA (punto extra)

### Uso
```bash
python scripts/ex4.py --input results/blast_results.txt --pattern "Homo sapiens" --output results/filtered_hits.txt [--fasta results/sequences.fasta] [--email email@example.com]
```

### Parámetros
- `--input/-i`: archivo de reporte BLAST (formato txt) - **REQUERIDO**
- `--pattern/-p`: patrón a buscar en las descripciones (ej: "Homo sapiens") - **REQUERIDO**
- `--output/-o`: archivo de salida con lista de hits filtrados - **REQUERIDO**
- `--fasta/-f`: archivo FASTA opcional para guardar secuencias completas de GenBank (punto extra)
- `--email/-e`: email para identificación con NCBI Entrez (default: "lrebuffo@frba.utn.edu.ar")
- `--delay`: delay entre consultas a GenBank en segundos (default: 0.5)

### Ejemplo
```bash
python scripts/ex4.py -i results/blast_results.txt -p "Homo sapiens" -o results/filtered_hits.txt --fasta results/homo_sapiens_sequences.fasta
```

## Flujo de Trabajo Completo

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

4. **Realizar alineamiento múltiple de secuencias:**
   ```bash
   python scripts/ex3.py -i ../data/secuencias -o results/alignment.aln --verbose
   ```

5. **Analizar reporte BLAST con filtrado por patrón:**
   ```bash
   python scripts/ex4.py -i results/blast_results.txt -p "Homo sapiens" -o results/filtered_hits.txt --fasta results/homo_sapiens_sequences.fasta
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

6. **Alineamiento Múltiple Avanzado**:
   - Auto-detección de Clustal Omega
   - Limpieza automática de secuencias
   - Normalización de IDs de secuencias
   - Manejo robusto de archivos de entrada

## Notas Importantes

- **Conexión a Internet:** El Ejercicio 2 requiere conexión a internet para acceder a los servidores BLAST de NCBI
- **Rate Limiting:** Se incluye un delay entre búsquedas BLAST para no sobrecargar los servidores de NCBI
- **Clustal Omega:** El Ejercicio 3 requiere Clustal Omega instalado (incluido en el proyecto o en PATH del sistema)
- **Formato de Salida:** Los resultados BLAST incluyen información detallada de alineamientos, valores E, identidad y descripciones de proteínas
- **Manejo de Errores:** Todos los scripts incluyen manejo robusto de errores y mensajes informativos
- **Limpieza de Datos:** El Ejercicio 3 limpia automáticamente las secuencias, manteniendo solo nucleótidos válidos
- **GenBank Access:** El Ejercicio 4 requiere conexión a internet para acceder a GenBank cuando se usa la opción `--fasta`
- **NCBI Email:** El Ejercicio 4 requiere un email válido para identificarse con NCBI Entrez

## Dependencias

- `biopython>=1.83,<2`: Para manipulación de secuencias y BLAST
- `python>=3.9`: Versión mínima de Python requerida
- `clustal-omega`: Para alineamiento múltiple de secuencias (incluido en el proyecto)

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

### Ejercicio 3 (scripts/ex3.py) - Arquitectura de Alineamiento Múltiple

El script implementa un pipeline completo de alineamiento múltiple:

1. **`parse_args()`**: Configuración de parámetros de entrada y salida
2. **`clean_sequence_id(record_id)`**: Limpia y normaliza IDs de secuencias
3. **`clean_sequence_nucleotides(sequence)`**: Limpia secuencias manteniendo solo nucleótidos válidos
4. **`find_fasta_files(input_dir)`**: Localiza archivos FASTA en el directorio especificado
5. **`combine_fasta_files(input_files, output_file)`**: Combina múltiples archivos FASTA
   - Procesa cada archivo FASTA individualmente
   - Limpia secuencias e IDs
   - Genera archivo FASTA combinado limpio
6. **`find_clustal_omega_executable(clustalo_path)`**: Auto-detecta Clustal Omega
   - Busca en directorio del script
   - Busca en PATH del sistema
   - Permite especificación manual
7. **`run_clustal_omega(input_file, output_file, clustalo_exe, verbose)`**: Ejecuta alineamiento
   - Configura parámetros de Clustal Omega
   - Maneja salida y errores
   - Proporciona feedback detallado
8. **`main()`**: Orquesta el pipeline completo
   - Validación de entradas
   - Manejo de errores robusto
   - Mensajes informativos consistentes

**Flujo de datos:**
```
Directorio FASTA → Combinación → Limpieza → Clustal Omega → Alineamiento
```

### Ejercicio 4 (scripts/ex4.py) - Arquitectura de Filtrado BLAST

El script implementa un pipeline de filtrado y obtención de secuencias:

1. **`parse_args()`**: Configuración de parámetros de entrada y filtrado
2. **`parse_blast_report(file_path)`**: Parsea archivo de reporte BLAST en formato texto
   - Detecta inicio de resultados para cada query
   - Extrae información de cada hit (ID, Descripción, ACCESSION, etc.)
   - Maneja diferentes formatos de ACCESSION (emb|ACCESSION|, ref|ACCESSION|, etc.)
3. **`filter_hits_by_pattern(hits, pattern, case_sensitive)`**: Filtra hits por patrón
   - Búsqueda case-insensitive por defecto
   - Busca el patrón en la descripción de cada hit
4. **`extract_accessions(hits)`**: Extrae accessions únicos de hits filtrados
5. **`fetch_sequences_from_genbank(accessions, email, delay)`**: Obtiene secuencias de GenBank
   - Usa Bio.Entrez para consultar GenBank
   - Implementa delays para no sobrecargar NCBI
   - Maneja errores de conexión y acceso
6. **`format_filtered_hits(hits, pattern)`**: Formatea hits filtrados para reporte
7. **`write_fasta_file(file_path, sequences)`**: Escribe secuencias en formato FASTA
8. **`main()`**: Orquesta el pipeline completo

**Flujo de datos:**
```
Reporte BLAST → Parseo → Filtrado por Patrón → Extracción ACCESSION → GenBank → FASTA
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

### Ejercicio 3 - Alineamiento Múltiple (Formato FASTA)
```
>Homo_sapiens
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>Pan_troglodytes
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>Pongo_abelii
ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
```

### Ejercicio 4 - Reporte de Hits Filtrados y Secuencias FASTA

**Salida de hits filtrados:**
```
Análisis de Reporte BLAST - Filtrado por Patrón
============================================================
Patrón buscado: 'Homo sapiens'
Total de hits encontrados: 21

=== Query: NM_000138.5 ===
Hits encontrados: 6

Hit #1:
  ID de Proteína: NP_001393645.1
  ACCESSION: NP_001393645.1
  Descripción: fibrillin 1 [Homo sapiens]
  Longitud: 2871 aa
  Valor E: 0.00e+00
  Puntuación: 14956.0
  ...
```

**Salida FASTA (opcional):**
```
>NP_001393645.1 fibrillin 1 [Homo sapiens]
MRRGRLLEIALGFTVLLASYTSHGADANLEAGNVKETRASRAKRRGGGGHDALKGPNVCG
...
```

## Contribuciones

Este proyecto forma parte de un trabajo práctico de bioinformática. Para modificaciones o mejoras, contactar al autor del proyecto.


