import requests
import re
import logging
import csv
from Bio import SeqIO  # To properly handle FASTA files

logger = logging.getLogger()

url = "https://rest.uniprot.org/uniprotkb/search"

def extract_protein_name(protein_data):
    """
    Extrae el mejor nombre posible para la proteína:
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
    Extrae el primer número EC si está presente.
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
    Consulta UniProt con accession para extraer cofactores dentro de comments.
    """
    url_entry = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    try:
        response = requests.get(url_entry)
        response.raise_for_status()
        entry = response.json()
    except Exception as e:
        print(f"❌ Error al obtener entrada UniProt para {accession}: {e}")
        return "Error"

    cofactors = []

    comments = entry.get("comments")
    if comments is None:
        print(f"⚠️ 'comments' no existe para {accession}")
        return "None"
    if not isinstance(comments, list):
        print(f"⚠️ 'comments' no es una lista para {accession}. Es: {type(comments)}")
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
    Extrae dominios Pfam (id y EntryName) desde una entrada de UniProt.
    """
    url_entry = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    try:
        response = requests.get(url_entry)
        response.raise_for_status()
        entry = response.json()
    except Exception as e:
        print(f"❌ Error al obtener Pfam para {accession}: {e}")
        return "Error"

    pfam_domains = []

    for xref in entry.get("uniProtKBCrossReferences", []):
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
    Extrae el ID de AlphaFoldDB desde una entrada de UniProt.
    """
    url_entry = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    try:
        response = requests.get(url_entry)
        response.raise_for_status()
        entry = response.json()
    except Exception as e:
        print(f"❌ Error al obtener AlphaFoldDB para {accession}: {e}")
        return "Error"

    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref.get("database") == "AlphaFoldDB":
            return xref.get("id", "NoID")

    return "None"

def annotate_uniprot_codes(
    valid_wp_codes,
    output_file="annotation.csv",
    incluir_organismo=True,
    incluir_nombre=True,
    incluir_ec=True,
    incluir_cofactores=True,
    incluir_pfam=True,
    incluir_alphafold=True
):
    results = []
    url = "https://rest.uniprot.org/uniprotkb/search"

    for wp in valid_wp_codes:
        query_string = f"xref:RefSeq-{wp}"
        params = {
            "fields": "accession,organism_name,protein_name",
            "query": query_string,
            "format": "json"
        }

        try:
            response = requests.get(url, params=params)
            response.raise_for_status()
            data = response.json()

            if data.get("results"):
                for r in data["results"]:
                    row = {"WP_code": wp}

                    accession = r.get("primaryAccession", "Not found")
                    row["UniProt_accession"] = accession

                    # Organismo
                    if incluir_organismo:
                        row["Organism"] = r.get("organism", {}).get("scientificName", "Not found")

                    # Nombre de proteína y EC
                    protein_data = r.get("proteinDescription", {})
                    if incluir_nombre:
                        row["Protein_name"] = extract_protein_name(protein_data)
                    if incluir_ec:
                        row["EC_number"] = extract_ec_number(protein_data)

                    # Datos adicionales con consulta por accession
                    if accession != "Not found":
                        if incluir_cofactores:
                            row["Cofactors"] = get_cofactors_from_accession(accession)
                        if incluir_pfam:
                            row["Pfam_domains"] = get_pfam_domains_from_accession(accession)
                        if incluir_alphafold:
                            row["AlphaFoldDB_ID"] = get_alphafold_id_from_accession(accession)
                    else:
                        if incluir_cofactores: row["Cofactors"] = "None"
                        if incluir_pfam: row["Pfam_domains"] = "None"
                        if incluir_alphafold: row["AlphaFoldDB_ID"] = "None"

                    results.append(row)
            else:
                # Sin resultados
                row = {"WP_code": wp, "UniProt_accession": "Not found"}
                if incluir_organismo: row["Organism"] = "Not found"
                if incluir_nombre: row["Protein_name"] = "Not found"
                if incluir_ec: row["EC_number"] = "None"
                if incluir_cofactores: row["Cofactors"] = "None"
                if incluir_pfam: row["Pfam_domains"] = "None"
                if incluir_alphafold: row["AlphaFoldDB_ID"] = "None"
                results.append(row)

        except Exception as e:
            print(f"❌ Error al consultar {wp}: {e}")
            row = {"WP_code": wp, "UniProt_accession": "error"}
            if incluir_organismo: row["Organism"] = "error"
            if incluir_nombre: row["Protein_name"] = "error"
            if incluir_ec: row["EC_number"] = "error"
            if incluir_cofactores: row["Cofactors"] = "error"
            if incluir_pfam: row["Pfam_domains"] = "error"
            if incluir_alphafold: row["AlphaFoldDB_ID"] = "error"
            results.append(row)

    # Determinar columnas a escribir dinámicamente
    all_possible_fields = [
        "WP_code", "UniProt_accession", "Organism", "Protein_name",
        "EC_number", "Cofactors", "Pfam_domains", "AlphaFoldDB_ID"
    ]
    included_fields = ["WP_code", "UniProt_accession"]
    if incluir_organismo: included_fields.append("Organism")
    if incluir_nombre: included_fields.append("Protein_name")
    if incluir_ec: included_fields.append("EC_number")
    if incluir_cofactores: included_fields.append("Cofactors")
    if incluir_pfam: included_fields.append("Pfam_domains")
    if incluir_alphafold: included_fields.append("AlphaFoldDB_ID")

    with open(output_file, mode="w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=included_fields, delimiter=';')
        writer.writeheader()
        writer.writerows(results)

    print(f"✅ Archivo CSV generado en: {output_file}")
