from Bio import SeqIO, AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Seq import Seq
import os, glob

print('BIOPYTHON VERSION :', __import__('Bio').__version__)

# Buscar archivos FASTA en el directorio data
input_files = glob.glob("../data/secuencias/*.fasta")
combined_fasta = "secuencias_combinadas.fa"
aln_file = "secuencias_alineadas.aln"

# Borrar archivos previos si existen
for f in [combined_fasta, aln_file]:
    if os.path.exists(f):
        os.remove(f)

# Combinar las secuencias en un solo FASTA limpio
with open(combined_fasta, "w", encoding="ascii", newline="\n") as output_handle:
    for file in input_files:
        for record in SeqIO.parse(file, "fasta"):
            clean_id = record.id.strip().replace(" ", "_")
            clean_seq = ''.join(ch for ch in str(record.seq).upper() if ch in "ACGTU")
            record.id = clean_id
            record.description = ""
            record.seq = Seq(clean_seq)
            SeqIO.write(record, output_handle, "fasta")

print(f"Archivo combinado creado: {combined_fasta}")
print("--------------")

script_dir = os.path.dirname(os.path.abspath(__file__))  # carpeta donde está el ejecutable de clustal omega
clustalo_exe = os.path.join(script_dir, "clustal-omega-1.2.2", "clustalo.exe")

# Configurar el comando
clustalo_cline = ClustalOmegaCommandline(
    cmd=clustalo_exe,
    infile=combined_fasta,
    outfile=aln_file,
    verbose=True,
    auto=True,
    force=True
)

#print("Comando Clustal Omega:", clustalo_cline)
#print("--------------")

try:
    # Ejecutar Clustal Omega
    stdout, stderr = clustalo_cline()
    print(stdout)
    print(stderr)
    print("Clustal Omega ejecutado correctamente")
    print("--------------")

    # Leer el alineamiento resultante
    if not os.path.exists(aln_file):
        raise FileNotFoundError(f"No se generó el archivo {aln_file}")

    #align = AlignIO.read(aln_file, "fasta")
    #print(align)
    #print("--------------")

    #for i, record in enumerate(align):
    #    print(f"{i+1} ---> {record.id} : {record.seq}")

except Exception as E:
    print("Error:", E)
    print("--------------")
