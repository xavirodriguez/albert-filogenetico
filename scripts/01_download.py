"""
HIV-1 Sequence Downloader from NCBI Entrez.

This module provides functionality to download HIV-1 sequences from the NCBI
Nucleotide database using specified search queries. It handles batch downloads
and saves the results in FASTA format.

Biological Context:
    The acquisition of high-quality genomic data is the foundation of any
    evolutionary analysis. For HIV-1, this typically involves retrieving
    complete genomes or specific genes (e.g., pol, env) from global repositories
    to capture the virus's immense genetic diversity.

Pipeline Stage:
    Phase 1 of 7: Data Acquisition.

Example:
    >>> # Run from terminal:
    >>> # python scripts/01_download.py
"""

from __future__ import annotations
import os
import yaml
import logging
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from datetime import datetime
from typing import Any, Dict, List

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("logs/download.log"),
        logging.StreamHandler()
    ]
)

def load_config() -> Dict[str, Any]:
    """
    Load pipeline configuration from a YAML file.

    Returns:
        Dict[str, Any]: Configuration dictionary containing paths and parameters.

    Raises:
        FileNotFoundError: If 'config.yaml' is missing.
    """
    with open("config.yaml", "r") as f:
        return yaml.safe_load(f)

def download_sequences(
    email: str,
    virus_query: str,
    max_sequences: int,
    output_file: str
) -> List[SeqRecord]:
    """
    Download sequences from NCBI Entrez and save them to a FASTA file.

    Args:
        email (str): Valid email address required by NCBI for Entrez usage.
        virus_query (str): Entrez search query string.
        max_sequences (int): Maximum number of sequences to retrieve.
        output_file (str): Path where the resulting FASTA file will be saved.

    Returns:
        List[SeqRecord]: A list of Biopython SeqRecord objects.

    Notes:
        NCBI rate limits apply. This function uses batching to handle large
        requests (100 sequences per batch).
    """
    Entrez.email = email
    logging.info(f"Searching for {virus_query} in NCBI...")
    
    handle = Entrez.esearch(db="nucleotide", term=virus_query, retmax=max_sequences)
    record = Entrez.read(handle)
    handle.close()
    
    id_list = record["IdList"]
    logging.info(f"Found {len(id_list)} sequences. Starting download...")
    
    sequences: List[SeqRecord] = []
    batch_size = 100
    for i in range(0, len(id_list), batch_size):
        batch_ids = id_list[i:i+batch_size]
        logging.info(f"Downloading batch {i//batch_size + 1}...")
        
        handle = Entrez.efetch(db="nucleotide", id=batch_ids, rettype="fasta", retmode="text")
        batch_records = list(SeqIO.parse(handle, "fasta"))
        sequences.extend(batch_records)
        handle.close()
    
    with open(output_file, "w") as f:
        SeqIO.write(sequences, f, "fasta")
    
    logging.info(f"Downloaded {len(sequences)} sequences to {output_file}")
    return sequences

def main() -> None:
    """
    Entry point for the download script.
    """
    config = load_config()
    raw_data_path = config["paths"]["raw_data"]
    os.makedirs(raw_data_path, exist_ok=True)
    
    # Query for HIV-1 complete genomes or large segments
    query = "HIV-1[Organism] AND (complete genome[Title] OR pol[Gene])"
    email = "your_email@example.com"
    
    output_file = os.path.join(raw_data_path, "hiv_raw.fasta")
    
    download_sequences(email, query, 100, output_file)

if __name__ == "__main__":
    main()
