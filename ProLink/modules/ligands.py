
import csv
import logging
import re
import requests
from Bio import SeqIO


logger = logging.getLogger()

def extract_pdb_codes_from_fasta(fasta_file:str) -> set:
    '''
    Extract PDB codes from the FASTA file
    
    Parameters
    ----------
    fasta_file : str
        Path to the input FASTA file
    
    Returns
    -------
    set
        A set of unique PDB codes found in the FASTA file
    '''
    logger.info("Attempting to extract PDB codes from FASTA")
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    pdb_codes = set()

    for seq in sequences:
        match = re.search(r'\b([A-Za-z0-9]{4})_[A-Za-z]\b', seq.description)
        if match:
            code = match.group(1).split("_")[0]  # Only the 4 characters before the "_"
            pdb_codes.add(code.upper())

    return pdb_codes

def get_ligands_from_pdb(pdb_code:str) -> list:
    '''
    Query the PDB API to extract ligands for a given code
    
    Parameters
    ----------
    pdb_code : str
        The PDB code to query
    
    Returns
    -------
    list of tuples
        A list of tuples containing (ligand ID, ligand name)
    '''
    base_url = "https://data.rcsb.org/rest/v1/core"
    entry_url = f"{base_url}/entry/{pdb_code}"

    response = requests.get(entry_url)
    if not response.ok:
        logger.error(f"ERROR: Could not access entry {pdb_code}")
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

def annotate_ligands_from_fasta(fasta_file:str, output_csv:str) -> None:
    '''
    Main function that extracts PDB codes and annotates their ligands into a CSV file
    
    Parameters
    ----------
    fasta_file : str
        Path to the input FASTA file
    output_csv : str
        Path to the output CSV file
    '''
    logger.info("Entering annotate_ligands_from_fasta")
    pdb_codes = extract_pdb_codes_from_fasta(fasta_file)
    logger.info(f"PDB codes found: {pdb_codes}")

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

    logger.debug(f"Output CSV path: {output_csv}")
    with open(output_csv, "w", newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile, delimiter=';')  # Separator compatible with Spanish Excel
        writer.writerow(headers)
        writer.writerows(all_data)

    logger.info(f"CSV file generated: {output_csv}")
