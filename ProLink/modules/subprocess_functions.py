
import logging
import subprocess
import time
import os
import re

from .. import ProLink_path

logger = logging.getLogger()

def clean_label(label, protein_name=""): 
    # Removes WP/XP/NP codes
    label = re.sub(r'(W|X|N)P[\s_]\d{9}\.\d', '', label)

    # Removes "MULTISPECIES:" and descriptions
    label = re.sub(r'MULTISPECIES:\s*', '', label, flags=re.IGNORECASE)

    # Removes the protein name if present
    if protein_name:
        protein_parts = re.split(r'[_\s]+', protein_name)
        protein_regex = r'[\s_\-]*'.join(map(re.escape, protein_parts))
        label = re.sub(protein_regex, "", label, flags=re.IGNORECASE)

    # Performs preliminary cleanup of specific words
    label = re.sub(r'unclassified', '', label, flags=re.IGNORECASE)
    label = re.sub(r'Same[\s_]+Domains', '', label, flags=re.IGNORECASE)

    # Removes isolated hyphens
    label = re.sub(r'[-]*', '', label).strip()

    # Removes leading and trailing quotation marks, if any
    label = label.strip("'\"")

    # Prevents issues with parentheses in Newick labels
    label = re.sub(r'[()]', '_', label)

    # Abbreviates the genus ONLY if there are at least two words and the second is not "sp."
    label = re.sub(
        r"^[\s_]*([A-Z])[a-zA-Z0-9]+[\s_]+(?!sp[\s\._])([a-z]+)",
        r"\1_\2",
        label
    )

    return label.strip(" _")


def clean_newick_string(newick_str, protein_name):
    """
    Cleans all labels in a Newick tree string by applying the clean_label function.
    """
    if not protein_name:
        raise ValueError("A protein name was expected but not received")
    print(f" [DEBUG] Protein name received clean_newick_string: {protein_name}")

    # Pattern to match labels that include the cluster marker. It matches labels that are
    # either quoted (single or double) or not quoted.
    pattern = re.compile(
        r"('([^']+---C\d+[^']*)'|\"([^\"]+---C\d+[^\"]*)\"|([A-Za-z0-9 _\.\-]+---C\d+))",
        flags=re.IGNORECASE
    )
    def replacer(match):
        full_label = match.group(0)
        return clean_label(full_label, protein_name)
    return pattern.sub(replacer, newick_str)

def align(muscle_input:str, muscle_output:str) -> None:
    '''
    Run a local alignment with MUSCLE v5
    Parameters
    ----------
    muscle_input : str
        Path of the input MUSCLE file
    muscle_output : str
        Path of the output MUSCLE file
    '''
    logging.info(f"\n-- Aligning sequences with MUSCLE")
    muscle_cmd = ['muscle', '-super5', muscle_input, '-output', muscle_output]
    logging.debug(f"Running MUSCLE alignment: {' '.join(muscle_cmd)}")
    muscle_run = subprocess.run(muscle_cmd)
    if muscle_run.returncode != 0:
        logger.error(f"ERROR: MUSCLE failed")
        raise RuntimeError(f"MUSCLE failed")

def tree(tree_type:str, bootstrap_replications:int, muscle_output:str, mega_output:str, protein_name:str) -> None:
    '''
    Run MEGA-CC to generate a phylogenetic tree
    Parameters
    ----------
    tree_type : str
        Type of tree to generate
    bootstrap_replications : int
        Number of bootstrap replications
    muscle_output : str
        Path of the input file (FASTA format from MUSCLE)
    mega_output : str
        Path of the MEGA-CC output file
    '''
    # Build the path to the MEGA-CC configuration file
    mega_config_input = f"{ProLink_path}/mega_configs/{tree_type}_{bootstrap_replications}.mao"
    logging.info(f"\n-- Generating phylogenetic tree with MEGA-CC")
    mega_cmd = ['megacc', '-a', mega_config_input, '-d', muscle_output, '-o', mega_output]
    logging.debug(f"Running MEGA-CC: {' '.join(mega_cmd)}")

    # Capture stdout and stderr to review MEGA-CC messages
    mega_run = subprocess.run(mega_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    logger.debug(f"MEGA-CC stdout: {mega_run.stdout}")
    logger.debug(f"MEGA-CC stderr: {mega_run.stderr}")
    if mega_run.returncode != 0:
        logger.error("ERROR: MEGA-CC failed")
        raise RuntimeError("MEGA-CC failed")

    # Wait for a few seconds to give time for the output file to be written
    time.sleep(5)

    # Verify that the output file exists before attempting to clean it
    if not os.path.exists(mega_output):
        # If the expected file with .nwk is not found, try the alternative with .mega
        alternative = mega_output.rsplit('.', 1)[0] + ".mega"
        if os.path.exists(alternative):
            logging.info(f"Using alternative output file: {alternative}")
            mega_output = alternative
        else:
            logger.error(f"ERROR: MEGA-CC did not produce the output file: {mega_output}")
            raise FileNotFoundError(f"Output file {mega_output} not found")

    # Read the generated Newick tree and clean its labels
    try:
        with open(mega_output, 'r') as f:
            newick = f.read()
        cleaned_newick = clean_newick_string(newick, protein_name=protein_name)
        with open(mega_output, 'w') as f:
            f.write(cleaned_newick)
        logging.info(f"Cleaned Newick tree saved in '{mega_output}'")
        logging.info("Tree cleaned and abbreviated")
    except Exception as e:
        logger.error(f"ERROR while cleaning the Newick file: {e}")
        raise
        
