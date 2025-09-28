## Ejercicio 1 – Procesamiento de secuencias (Python)

Este script lee uno o más archivos GenBank con mRNA(s) de referencia (NM_...), genera las 6 traducciones posibles (3 marcos forward y 3 reverse complement), detecta ORFs y escribe un FASTA con las secuencias de aminoácidos. Marca como "best_frame=true" al ORF más largo encontrado entre los 6 marcos.

### Requisitos
- Python 3.9+
- Dependencias: ver `requirements.txt`

Instalar dependencias:
```bash
pip install -r requirements.txt
```

### Uso
```bash
python ex1.py --input NMxxxx.gbk NMyyyy.gbk --output out.fasta [--min_orf_len 30] [--write_all_frames]
```

Parámetros:
- `--input/-i`: uno o más archivos GenBank (`.gb`/`.gbk`) con mRNA(s) NM_. Si un archivo contiene múltiples registros, se procesan todos.
- `--output/-o`: archivo FASTA de salida con aminoácidos.
- `--min_orf_len`: longitud mínima de ORF en aa (default 30).
- `--write_all_frames`: si se incluye, agrega también las traducciones completas de los 6 marcos.

### Notas
- El mejor marco se define como aquel con el ORF de mayor longitud. Si prefieres otra heurística (p. ej., presencia de CDS anotado), avísame y lo ajustamos.
- Para obtener el input desde NCBI: busca el gen (p. ej., INS), filtra por RefSeq mRNA (NM_), y exporta en formato GenBank.

### Ejemplo mínimo
```bash
python ex1.py -i NM_000207.gbk -o INS_orfs.fasta --min_orf_len 30
```


