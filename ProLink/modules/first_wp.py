
import logging
import requests
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


logger = logging.getLogger()

def get_wp_from_code(code: str) -> str:
    """
    Get the WP code from an input code
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"({code})",
        "fields": "accession",
        "format": "json"
    }

    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        data = response.json()

        if not data.get("results"):
            return None

        accession = data["results"][0].get("primaryAccession", None)
        if not accession:
            return None

        # Get the WP from the input
        entry_url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
        entry_resp = requests.get(entry_url)
        entry_resp.raise_for_status()
        entry = entry_resp.json()

        for ref in entry.get("uniProtKBCrossReferences", []):
            if ref.get("database") == "RefSeq":
                wp_of_interest = ref.get("id")
                if wp_of_interest and wp_of_interest.startswith("WP_"):
                    return wp_of_interest
        return None

    except Exception as e:
        logger.error(f"ERROR: Error retrieving WP from {code}: {e}")
        return None

def reorder_fasta_with_study_sequence(txt_file, fasta_file, study_wp, study_fasta_file, output_fasta_file):
    # Read the txt file to find the study WP cluster and its descriptions
    cluster_dict = {}
    current_cluster = None
    study_header = None

    with open(txt_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#Cluster'):
                current_cluster = line.replace('#Cluster', '').strip()
                cluster_dict[current_cluster] = []
            elif line:
                parts = line.split(maxsplit=1)
                wp = parts[0]
                desc = parts[1] if len(parts) > 1 else ""
                cluster_dict[current_cluster].append((wp, desc))
                if wp == study_wp:
                    study_header = f"{wp} {desc}---C{current_cluster}"
                    study_cluster = current_cluster

    if study_header is None:
        logger.error(f"ERROR: Code {study_wp} not found in the txt file")
        return

    # Read the study sequence from my_sequence.fasta
    study_seq_record = next(SeqIO.parse(study_fasta_file, "fasta"))
    study_seq = str(study_seq_record.seq)
    new_study_record = SeqRecord(
        seq=study_seq_record.seq,
        id="",  # Remove automatic ID
        description=study_header
    )

    # Read the original fasta file and sort sequences by cluster
    fasta_records_by_cluster = {cluster: [] for cluster in cluster_dict}
    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.description
        if '---C' in header:
            cluster_id = header.split('---C')[-1]
            if cluster_id in fasta_records_by_cluster:
                fasta_records_by_cluster[cluster_id].append(record)

    # Replace the first sequence of the corresponding cluster
    original_records = fasta_records_by_cluster[study_cluster]
    updated_records = [new_study_record] + original_records[1:]  # Replace header
    fasta_records_by_cluster[study_cluster] = updated_records

    # Write all sequences to the new fasta file
    all_new_records = []
    for cluster_id in sorted(fasta_records_by_cluster.keys(), key=lambda x: int(x)):
        all_new_records.extend(fasta_records_by_cluster[cluster_id])

    SeqIO.write(all_new_records, output_fasta_file, "fasta")
    logger.info(f"Updated fasta file: {output_fasta_file}")
