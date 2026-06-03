"""
Quality Control and Filtering for HIV-1 Sequences.

This module performs biological quality control on raw FASTA sequences,
filtering out entries that do not meet length requirements or contain
excessive ambiguous bases (Ns).

Biological Context:
    In HIV-1 evolutionary studies, sequences with high ambiguity or short
    lengths can introduce noise in phylogenetic reconstruction. Standard
    practice involves filtering sequences that cover less than 80% of the
    target region or have >5% ambiguous sites.

Pipeline Stage:
    Phase 2 of 7: Quality Control.

Example:
    >>> # Run from terminal:
    >>> # python scripts/02_qc.py
"""

from __future__ import annotations
import os
import yaml
import logging
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from typing import Any, Dict

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("logs/qc.log"),
        logging.StreamHandler()
    ]
)

def load_config() -> Dict[str, Any]:
    """
    Load pipeline configuration from a YAML file.

    Returns:
        Dict[str, Any]: Configuration dictionary.
    """
    with open("config.yaml", "r") as f:
        return yaml.safe_load(f)

def run_qc(
    input_file: str,
    min_len: int,
    max_ns: float,
    processed_path: str
) -> pd.DataFrame:
    """
    Execute biological quality control filtering on a FASTA file.

    Args:
        input_file (str): Path to the input raw FASTA file.
        min_len (int): Minimum sequence length required (in bp).
        max_ns (float): Maximum allowed percentage of ambiguous bases (Ns).
        processed_path (str): Directory where filtered output and stats will be saved.

    Returns:
        pd.DataFrame: Summary statistics of all analyzed sequences.

    Notes:
        The function saves `hiv_filtered.fasta` and `qc_stats.csv` to the
        `processed_path` directory.
    """
    logging.info(f"Starting Quality Control on {input_file}...")
    
    stats = []
    filtered_sequences = []
    
    for record in SeqIO.parse(input_file, "fasta"):
        seq_len = len(record.seq)
        n_count = record.seq.upper().count('N')
        ns_percent = (n_count / seq_len) * 100
        gc = gc_fraction(record.seq) * 100
        
        stats.append({
            'id': record.id,
            'length': seq_len,
            'ns_percent': ns_percent,
            'gc_content': gc
        })
        
        # Rigorous biological filtering
        if seq_len >= min_len and ns_percent <= max_ns:
            filtered_sequences.append(record)
            
    df = pd.DataFrame(stats)
    
    # Save statistics
    stats_file = os.path.join(processed_path, "qc_stats.csv")
    df.to_csv(stats_file, index=False)
    
    # Save filtered sequences
    output_fasta = os.path.join(processed_path, "hiv_filtered.fasta")
    with open(output_fasta, "w") as f:
        SeqIO.write(filtered_sequences, f, "fasta")
        
    logging.info(f"QC Complete. {len(filtered_sequences)} sequences passed filters.")
    
    return df

def generate_plots(df: pd.DataFrame, output_path: str) -> None:
    """
    Generate distribution plots for sequence metrics.

    Args:
        df (pd.DataFrame): DataFrame containing 'length' and 'gc_content' columns.
        output_path (str): Directory where the PNG plot will be saved.
    """
    plt.figure(figsize=(12, 6))
    
    plt.subplot(1, 2, 1)
    sns.histplot(df['length'], kde=True, color='blue')
    plt.title('Distribution of Sequence Lengths')
    
    plt.subplot(1, 2, 2)
    sns.histplot(df['gc_content'], kde=True, color='green')
    plt.title('GC Content Distribution')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, "qc_metrics.png"))
    logging.info(f"QC plots saved to {output_path}")

def main() -> None:
    """
    Entry point for the QC script.
    """
    config = load_config()
    raw_file = os.path.join(config["paths"]["raw_data"], "hiv_raw.fasta")
    processed_path = config["paths"]["processed_data"]
    figures_path = config["paths"]["figures"]
    
    os.makedirs(processed_path, exist_ok=True)
    os.makedirs(figures_path, exist_ok=True)
    
    if not os.path.exists(raw_file):
        logging.error(f"Input file {raw_file} not found. Run download script first.")
        return

    df = run_qc(raw_file, config["min_length"], config["max_ns_percent"], processed_path)
    generate_plots(df, figures_path)

if __name__ == "__main__":
    main()
