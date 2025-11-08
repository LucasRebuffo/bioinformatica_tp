import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
input_fasta = os.path.join(script_dir, "..", "data", "NM_000138.5.fasta")
orfs_fasta = os.path.join(script_dir, "..", "results", "orfs_NM_000138.5.fasta")
proteins_fasta = os.path.join(script_dir, "..", "results", "proteins_NM_000138.5.fasta")
output_file = os.path.join(script_dir, "..", "results", "domains_NM_000138.5.txt")

# ----------------------------
# Ejecutar getorf de EMBOSS con subprocess
# ----------------------------

# Comando que se ejecutaría en la terminal
cmd = [
    "getorf",
    "-sequence", input_fasta,
    "-outseq", orfs_fasta,
    "-find", "3"  # buscar ORFs en ambas hebras
]

# Ejecutar el comando
try:
    subprocess.run(cmd, check=True)
    print("Ejecucion de EMBOSS GETORF completada.\n")
except FileNotFoundError:
    print("No se encontró el comando 'getorf'. Asegúrate de que EMBOSS esté instalado y en el PATH.")
    exit()
except subprocess.CalledProcessError as e:
    print(f"Error ejecutando getorf: {e}")
    exit()

# Verificar existencia de archivo
for f in [orfs_fasta]:
    if not os.path.exists(f):
        print(f"No se encontró el archivo: {os.path.abspath(f)}")
        exit()
    else:
        print(f"Archivo encontrado: {os.path.abspath(f)}")

# ----------------------------
# Traducir los ORFs a proteínas (necesario para patmatmotifs)
# ----------------------------

records = []
for record in SeqIO.parse(orfs_fasta, "fasta"):
    protein_seq = record.seq.translate(to_stop=True)
    new_record = SeqRecord(
        protein_seq,
        id=record.id,
        description="Translated protein from ORF"
    )
    records.append(new_record)

SeqIO.write(records, proteins_fasta, "fasta")
print(f"Archivo de proteinas traducidas generado: {proteins_fasta}")

# ----------------------------
# Ejecutar patmatmotifs (EMBOSS)
# ----------------------------
cmd = [
    "patmatmotifs",
    "-sequence", proteins_fasta,
    "-full",
    "-outfile", output_file
]

try:
    subprocess.run(cmd, check=True)
    print("Analisis de dominios completado con exito.\n")
except FileNotFoundError:
    print("No se encontró el comando 'patmatmotifs'. Asegúrate de que EMBOSS esté instalado y en el PATH.")
    exit()
except subprocess.CalledProcessError as e:
    print(f"Error ejecutando patmatmotifs: {e}")
    exit()
