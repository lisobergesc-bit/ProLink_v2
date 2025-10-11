
import csv
import logging
import requests
from Bio import SeqIO  # To properly handle FASTA files


logger = logging.getLogger()

url = "https://rest.uniprot.org/uniprotkb/search"

def extract_protein_name(protein_data):
    """
    Extracts the best possible protein name:
    1. recommendedName
    2. submissionNames[0]
    3. alternativeNames[0]
    """
    try:
        return protein_data["recommendedName"]["fullName"]["value"]
    except (KeyError, TypeError):
        pass
    try:
        return protein_data["submissionNames"][0]["fullName"]["value"]
    except (KeyError, IndexError, TypeError):
        pass
    try:
        return protein_data["alternativeNames"][0]["fullName"]["value"]
    except (KeyError, IndexError, TypeError):
        pass
    return "Not found"

def extract_ec_number(protein_data):
    """
    Extracts the first EC number if present
    """
    try:
        ec_list = protein_data.get("recommendedName", {}).get("ecNumbers", [])
        if ec_list:
            return ec_list[0]["value"]
    except (KeyError, IndexError, TypeError):
        pass
    return "None"

def get_cofactors_from_accession(accession):
    """
    Queries UniProt using the accession to extract cofactors from 'comments'
    """
    url_entry = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    try:
        response = requests.get(url_entry)
        response.raise_for_status()
        entry = response.json()
        logger.debug(f"UniProt entry (cofactors) for {accession} successfully retrieved")
    except Exception as e:
        logger.error(f"ERROR: Error retrieving UniProt entry for {accession}: {e}")
        return "Error"

    cofactors = []

    comments = entry.get("comments")
    logger.debug(f"Type of 'comments' for {accession}: {type(comments)}")

    if comments is None:
        logger.warning(f"WARNING: 'comments' does not exist for {accession}")
        return "None"
    if not isinstance(comments, list):
        logger.warning(f"WARNING: 'comments' is not a list for {accession}. Type: {type(comments)}")
        return "None"

    for comment in comments:
        if comment.get("commentType") == "COFACTOR":
            for cofactor in comment.get("cofactors", []):
                name_field = cofactor.get("name")
                if isinstance(name_field, dict):
                    name = name_field.get("value")
                elif isinstance(name_field, str):
                    name = name_field
                else:
                    name = None
                if name:
                    cofactors.append(name)

    return "; ".join(cofactors) if cofactors else "None"

def get_pfam_domains_from_accession(accession):
    """
    Extracts Pfam domains (id and EntryName) from a UniProt entry
    """
    url_entry = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    try:
        response = requests.get(url_entry)
        response.raise_for_status()
        entry = response.json()
        logger.debug(f"UniProt entry (Pfam) for {accession} successfully retrieved")
    except Exception as e:
        logger.error(f"ERROR: Error retrieving Pfam for {accession}: {e}")
        return "Error"

    pfam_domains = []

    crossrefs = entry.get("uniProtKBCrossReferences", [])
    logger.debug(f"Type of 'uniProtKBCrossReferences' for {accession}: {type(crossrefs)}")

    for xref in crossrefs:
        if xref.get("database") == "Pfam":
            pfam_id = xref.get("id", "NoID")
            entry_name = None
            for prop in xref.get("properties", []):
                if prop.get("key") == "EntryName":
                    entry_name = prop.get("value", "NoEntryName")
                    break
            pfam_domains.append(f"{pfam_id} ({entry_name})" if entry_name else pfam_id)

    return "; ".join(pfam_domains) if pfam_domains else "None"

def get_alphafold_id_from_accession(accession):
    """
    Extracts AlphaFoldDB ID from a UniProt entry
    """
    url_entry = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    try:
        response = requests.get(url_entry)
        response.raise_for_status()
        entry = response.json()
        logger.debug(f"UniProt entry (AlphaFold) for {accession} successfully retrieved")
    except Exception as e:
        logger.error(f"ERROR: Error retrieving AlphaFoldDB for {accession}: {e}")
        return "Error"

    crossrefs = entry.get("uniProtKBCrossReferences", [])
    logger.debug(f"Type of 'uniProtKBCrossReferences' for AlphaFold {accession}: {type(crossrefs)}")

    for xref in crossrefs:
        if xref.get("database") == "AlphaFoldDB":
            return xref.get("id", "NoID")

    return "None"

def annotate_uniprot_codes(
    valid_wp_codes,
    output_file="annotation.csv",
    include_organism=True,
    include_name=True,
    include_ec=True,
    include_cofactors=True,
    include_pfam=True,
    include_alphafold=True
):
    results = []
    url = "https://rest.uniprot.org/uniprotkb/search"

    logger.info(f"Starting annotation of {len(valid_wp_codes)} WP codes")

    for wp in valid_wp_codes:
        logger.info(f"\nQuerying UniProt for: {wp}")
        query_string = f"xref:RefSeq-{wp}"
        params = {
            "fields": "accession,organism_name,protein_name",
            "query": query_string,
            "format": "json"
        }

        try:
            response = requests.get(url, params=params)
            logger.debug(f"Queried URL: {response.url}")
            response.raise_for_status()
            data = response.json()
            logger.debug(f"JSON response for {wp}: {data}")
        except Exception as e:
            logger.error(f"ERROR: Error querying {wp}: {e}")
            row = {"WP_code": wp, "UniProt_accession": "error"}
            if include_organism: row["Organism"] = "error"
            if include_name: row["Protein_name"] = "error"
            if include_ec: row["EC_number"] = "error"
            if include_cofactors: row["Cofactors"] = "error"
            if include_pfam: row["Pfam_domains"] = "error"
            if include_alphafold: row["AlphaFoldDB_ID"] = "error"
            results.append(row)
            continue

        if data.get("results"):
            for r in data["results"]:
                logger.debug(f"UniProt entry processed for {wp}: {r}")
                row = {"WP_code": wp}
                accession = r.get("primaryAccession", "Not found")
                row["UniProt_accession"] = accession

                if include_organism:
                    row["Organism"] = r.get("organism", {}).get("scientificName", "Not found")

                protein_data = r.get("proteinDescription", {})
                logger.debug(f"proteinDescription for {wp}: {protein_data}")

                if include_name:
                    logger.info(f"Searching for protein name")
                    row["Protein_name"] = extract_protein_name(protein_data)
                if include_ec:
                    logger.info(f"Searching for EC number")
                    row["EC_number"] = extract_ec_number(protein_data)

                if accession != "Not found":
                    if include_cofactors:
                        logger.info(f"Searching for cofactors")
                        row["Cofactors"] = get_cofactors_from_accession(accession)
                    if include_pfam:
                        logger.info(f"Searching for Pfam domains")
                        row["Pfam_domains"] = get_pfam_domains_from_accession(accession)
                    if include_alphafold:
                        logger.info(f"Searching for AlphaFoldDB ID")
                        row["AlphaFoldDB_ID"] = get_alphafold_id_from_accession(accession)
                else:
                    if include_cofactors: row["Cofactors"] = "None"
                    if include_pfam: row["Pfam_domains"] = "None"
                    if include_alphafold: row["AlphaFoldDB_ID"] = "None"

                results.append(row)
        else:
            logger.warning(f"No results found for {wp}")
            row = {"WP_code": wp, "UniProt_accession": "Not found"}
            if include_organism: row["Organism"] = "Not found"
            if include_name: row["Protein_name"] = "Not found"
            if include_ec: row["EC_number"] = "None"
            if include_cofactors: row["Cofactors"] = "None"
            if include_pfam: row["Pfam_domains"] = "None"
            if include_alphafold: row["AlphaFoldDB_ID"] = "None"
            results.append(row)

    all_possible_fields = [
        "WP_code", "UniProt_accession", "Organism", "Protein_name",
        "EC_number", "Cofactors", "Pfam_domains", "AlphaFoldDB_ID"
    ]
    included_fields = ["WP_code", "UniProt_accession"]
    if include_organism: included_fields.append("Organism")
    if include_name: included_fields.append("Protein_name")
    if include_ec: included_fields.append("EC_number")
    if include_cofactors: included_fields.append("Cofactors")
    if include_pfam: included_fields.append("Pfam_domains")
    if include_alphafold: included_fields.append("AlphaFoldDB_ID")

    with open(output_file, mode="w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=included_fields, delimiter=';')
        writer.writeheader()
        writer.writerows(results)

    logger.info(f"CSV file generated at: {output_file}")
