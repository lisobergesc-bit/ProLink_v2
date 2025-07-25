import logging
import re
import requests
import csv
from Bio import SeqIO

logger = logging.getLogger()

def extract_pdb_codes_from_fasta(fasta_file):
    """Extrae códigos PDB del archivo FASTA"""
    logger.info("Intentando extraer códigos PDB del FASTA")
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    pdb_codes = set()

    for seq in sequences:
        match = re.search(r'\b([A-Za-z0-9]{4})_[A-Za-z]\b', seq.description)
        if match:
            code = match.group(1).split("_")[0]  # Solo los 4 caracteres antes del "_"
            pdb_codes.add(code.upper())

    return pdb_codes


def get_ligands_from_pdb(pdb_code):
    """Consulta la API del PDB para extraer ligandos de un código"""
    base_url = "https://data.rcsb.org/rest/v1/core"
    entry_url = f"{base_url}/entry/{pdb_code}"

    response = requests.get(entry_url)
    if not response.ok:
        print(f"No se pudo acceder a la entrada {pdb_code}")
        return []

    data = response.json()
    ligand_ids = data.get("rcsb_entry_container_identifiers", {}).get("non_polymer_entity_ids", [])

    ligands = []
    for ligand_id in ligand_ids:
        ligand_url = f"{base_url}/nonpolymer_entity/{pdb_code}/{ligand_id}"
        ligand_response = requests.get(ligand_url)
        if not ligand_response.ok:
            continue

        ligand_data = ligand_response.json()
        comp_id = ligand_data.get("pdbx_entity_nonpoly", {}).get("comp_id", "")
        name = ligand_data.get("pdbx_entity_nonpoly", {}).get("name", "")
        ligands.append((comp_id, name))

    return ligands


def annotate_ligands_from_fasta(fasta_file, output_csv):
    """Función final que extrae códigos PDB y anota sus ligandos en un CSV"""
    logger.info("Entrando en annotate_ligands_from_fasta")
    pdb_codes = extract_pdb_codes_from_fasta(fasta_file)
    print(f"Códigos PDB encontrados: {pdb_codes}")

    all_data = []
    max_ligands = 0

    for pdb in pdb_codes:
        ligands = get_ligands_from_pdb(pdb)
        row = [pdb]
        for comp_id, name in ligands:
            row.extend([comp_id, name])
        all_data.append(row)
        max_ligands = max(max_ligands, len(ligands))

    headers = ["PDB code"]
    for i in range(1, max_ligands + 1):
        headers.extend([f"Ligand {i}", f"Name {i}"])

    logger.debug(f"Ruta del CSV de salida: {output_csv}")
    with open(output_csv, "w", newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile, delimiter=';')  # Separador compatible con Excel español
        writer.writerow(headers)
        writer.writerows(all_data)

    print(f"✅ Archivo CSV generado: {output_csv}")
