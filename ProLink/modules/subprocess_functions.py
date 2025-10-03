
import logging
import subprocess
import time
import os
import re

from .. import ProLink_path

logger = logging.getLogger()

def clean_label(label, protein_name=""): 
    # Elimina códigos WP/XP/NP
    label = re.sub(r'(W|X|N)P[\s_]\d{9}\.\d', '', label)

    # Elimina "MULTISPECIES:" y descripciones
    label = re.sub(r'MULTISPECIES:\s*', '', label, flags=re.IGNORECASE)

    # Elimina nombre de la proteína si está presente
    if protein_name:
        protein_parts = re.split(r'[_\s]+', protein_name)
        protein_regex = r'[\s_\-]*'.join(map(re.escape, protein_parts))
        label = re.sub(protein_regex, "", label, flags=re.IGNORECASE)

    # Limpieza previa de palabras específicas
    label = re.sub(r'unclassified', '', label, flags=re.IGNORECASE)
    label = re.sub(r'Same[\s_]+Domains', '', label, flags=re.IGNORECASE)

    # Elimina guiones sueltos
    label = re.sub(r'[-]*', '', label).strip()

    # Elimina comillas iniciales y finales si existen
    label = label.strip("'\"")

    # ⚠️ Evita problemas de paréntesis en etiquetas Newick
    label = re.sub(r'[()]', '_', label)

    # Abrevia el género SOLO si hay al menos dos palabras y la segunda no es "sp."
    label = re.sub(
        r"^[\s_]*([A-Z])[a-zA-Z0-9]+[\s_]+(?!sp[\s\._])([a-z]+)",
        r"\1_\2",
        label
    )

    return label.strip(" _")


def clean_newick_string(newick_str, protein_name):
    if not protein_name:
        raise ValueError("❌ Se esperaba un nombre de proteína pero no ha llegado.")
    print(f" [DEBUG] Nombre de_la proteína recibido clean_newick_string: {protein_name}")
    pattern = re.compile(
        r"('([^']+---C\d+[^']*)'|\"([^\"]+---C\d+[^\"]*)\"|([A-Za-z0-9 _\.\-]+---C\d+))",
        flags=re.IGNORECASE
    )
    def replacer(match):
        full_label = match.group(0)
        return clean_label(full_label, protein_name)
    return pattern.sub(replacer, newick_str)

def align(muscle_input:str, muscle_output:str) -> None:
    logging.info(f"\n-- Aligning sequences with MUSCLE")
    muscle_cmd = ['muscle', '-super5', muscle_input, '-output', muscle_output]
    logging.debug(f"Running MUSCLE alignment: {' '.join(muscle_cmd)}")
    muscle_run = subprocess.run(muscle_cmd)
    if muscle_run.returncode != 0:
        logger.error(f"ERROR: MUSCLE failed")
        raise RuntimeError(f"MUSCLE failed")

def tree(tree_type:str, bootstrap_replications:int, muscle_output:str, mega_output:str, protein_name:str) -> None:
    mega_config_input = f"{ProLink_path}/mega_configs/{tree_type}_{bootstrap_replications}.mao"
    logging.info(f"\n-- Generating phylogenetic tree with MEGA-CC")
    mega_cmd = ['megacc', '-a', mega_config_input, '-d', muscle_output, '-o', mega_output]
    logging.debug(f"Running MEGA-CC: {' '.join(mega_cmd)}")

    mega_run = subprocess.run(mega_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    logger.debug(f"MEGA-CC stdout: {mega_run.stdout}")
    logger.debug(f"MEGA-CC stderr: {mega_run.stderr}")
    if mega_run.returncode != 0:
        logger.error("ERROR: MEGA-CC failed")
        raise RuntimeError("MEGA-CC failed")

    time.sleep(5)

    if not os.path.exists(mega_output):
        alternative = mega_output.rsplit('.', 1)[0] + ".mega"
        if os.path.exists(alternative):
            logging.info(f"Using alternative output file: {alternative}")
            mega_output = alternative
        else:
            logger.error(f"ERROR: MEGA-CC did not produce the output file: {mega_output}")
            raise FileNotFoundError(f"Output file {mega_output} not found")

    try:
        with open(mega_output, 'r') as f:
            newick = f.read()
        cleaned_newick = clean_newick_string(newick, protein_name=protein_name)
        with open(mega_output, 'w') as f:
            f.write(cleaned_newick)
        logging.info(f"Cleaned Newick tree saved in '{mega_output}'")
        logging.info("cleaned and abb")
    except Exception as e:
        logger.error(f"ERROR while cleaning the Newick file: {e}")
        raise
        
