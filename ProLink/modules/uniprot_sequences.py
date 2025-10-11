
import logging
import re
import requests
from Bio import SeqIO  # To properly handle FASTA files


logger = logging.getLogger()

url = "https://rest.uniprot.org/uniprotkb/search"

def check_uniprot_single(wp_code):
    """
    Verify the existence of a single WP code in UniProt.

    Parameters:
    wp_code (str): WP code to verify.

    Returns:
    bool: True if the WP code exists in UniProt, False otherwise.
    """
    params = {
        "query": f"xref:RefSeq-{wp_code}",
        "fields": "accession",
        "format": "json",
        "size": 1  # We only need to check if it exists
    }

    logger.debug(f"Checking to UniProt: {params['query']}")

    try:
        response = requests.get(url, params=params, timeout=10)
        response.raise_for_status()
        data = response.json()
        return bool(data.get("results"))  # Returns True if there's at least one result
    except requests.exceptions.RequestException as e:
        logger.error(f"ERROR: Could not connect to UniProt: {e}")
        return False

def filter_valid_sequences(input_fasta, output_fasta):
    """
    Filters sequences by removing those whose WP codes do not exist in UniProt.
    Sequences without a WP_ code are retained.

    Parameters:
    input_fasta (str): Input FASTA file with sequences.
    output_fasta (str): Output FASTA file with valid sequences.
    """
    sequences = list(SeqIO.parse(input_fasta, "fasta"))

    # Extract WP_ codes from sequence descriptions
    wp_data = {}
    for seq in sequences:
        match = re.search(r'((?:WP|XP|NP)_\d{9}\.\d)', seq.description)
        if match:
            wp_data[seq.description] = match.group(1)

    logger.debug(f"Extracted WP codes: {list(wp_data.values())}")

    logger.info(f"Total number of sequences: {len(sequences)}")
    logger.info(f"Number of WP codes found: {len(wp_data)}")

    # Verify each WP code in UniProt individually
    valid_wp_codes = {wp for wp in wp_data.values() if check_uniprot_single(wp)}

    # Filter valid sequences
    valid_sequences = [
        seq for seq in sequences
        if seq.description not in wp_data or wp_data[seq.description] in valid_wp_codes
    ]

    # Write the valid sequences to the new FASTA file
    SeqIO.write(valid_sequences, output_fasta, "fasta")
    logger.debug(f"Valid sequences after filtering: {len(valid_sequences)}")
    logger.info(f"Results saved to {output_fasta}")

    return valid_wp_codes
