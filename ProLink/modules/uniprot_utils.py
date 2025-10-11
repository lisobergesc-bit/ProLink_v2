
import logging
import requests
from xml.etree import ElementTree


logger = logging.getLogger()

def get_protein_name_from_wp(wp_code):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    
    # Get UID (NCBI ID)
    esearch_url = f"{base_url}esearch.fcgi?db=protein&term={wp_code}&retmode=xml"
    esearch_response = requests.get(esearch_url)
    esearch_root = ElementTree.fromstring(esearch_response.content)
    uid = esearch_root.findtext(".//Id")

    if not uid:
        logger.error(f"ERROR: No UID found for WP code {wp_code}")
        raise ValueError(f"No title found for UID {uid}")

    # Get name (Title)
    esummary_url = f"{base_url}esummary.fcgi?db=protein&id={uid}&retmode=xml"
    esummary_response = requests.get(esummary_url)
    esummary_root = ElementTree.fromstring(esummary_response.content)
    title = esummary_root.findtext(".//Item[@Name='Title']")

    if not title:
        logger.error(f"ERROR: No title found for UID {uid}")
        raise ValueError(f"No title found for UID {uid}")

    # Remove species name enclosed in brackets
    title = title.split(" [")[0].strip()

    # Remove only non-informative prefixes if they appear before a colon
    if ":" in title:
        prefix, rest = title.split(":", 1)
        prefix_clean = prefix.strip().upper()
        unwanted_prefixes = {
            "MULTISPECIES",
            "HYPOTHETICAL PROTEIN",
            "PUTATIVE",
            "UNNAMED PROTEIN PRODUCT",
            "UNCHARACTERIZED PROTEIN",
            "PREDICTED",
            "PROBABLE",
            "POSSIBLE",
            "GENERIC",
            "UNSPECIFIED",
            "UNKNOWN",
            "AUTOMATIC ANNOTATION",
            "PARTIAL",
            "LOW QUALITY PROTEIN"
        }
        if prefix_clean in unwanted_prefixes:
            logger.info(f"Prefix removed from title: '{prefix.strip()}'")
            title = rest.strip()

    return title
