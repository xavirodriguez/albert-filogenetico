"""
Multiple Sequence Alignment using MAFFT for HIV-1.

This module wraps the MAFFT alignment tool to align filtered HIV-1 sequences.
It is optimized for large datasets and maintains biological accuracy.

Biological Context:
    Accurate alignment is critical for phylogenetic inference. HIV-1 is
    characterized by hypervariable regions (especially in env) and frequent
    indels. MAFFT's FFT-NS-2 or PartTree algorithms are preferred for
    maintaining speed without sacrificing too much accuracy in large datasets.

Pipeline Stage:
    Phase 3 of 7: Sequence Alignment.

Example:
    >>> # Run from terminal:
    >>> # python scripts/03_align.py
"""

from __future__ import annotations
import os
import yaml
import subprocess
import logging
from typing import Any, Dict

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("logs/alignment.log"),
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

def run_mafft(input_file: str, output_file: str, config: Dict[str, Any]) -> None:
    """
    Execute MAFFT alignment with configured parameters.

    Args:
        input_file (str): Path to the input filtered FASTA file.
        output_file (str): Path where the aligned FASTA will be saved.
        config (Dict[str, Any]): Configuration dictionary containing tool paths.

    Notes:
        Requires MAFFT to be installed at the path specified in `config.yaml`.
        The `MAFFT_BINARIES` environment variable is automatically configured.
    """
    mafft_exe = config["alignment"]["executable"]
    mafft_binaries = config["alignment"]["binaries"]
    args = config["alignment"]["args"].split()
    
    logging.info(f"Running MAFFT alignment on {input_file}...")
    
    env = os.environ.copy()
    env["MAFFT_BINARIES"] = os.path.abspath(mafft_binaries)
    
    cmd = [mafft_exe] + args + [input_file]
    
    try:
        with open(output_file, "w") as out:
            result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True, env=env)
            
        if result.returncode != 0:
            logging.error(f"MAFFT error: {result.stderr}")
        else:
            logging.info(f"Alignment completed successfully. Saved to {output_file}")
            
    except Exception as e:
        logging.error(f"Failed to run MAFFT: {str(e)}")

def main() -> None:
    """
    Entry point for the alignment script.
    """
    config = load_config()
    input_file = os.path.join(config["paths"]["processed_data"], "hiv_filtered.fasta")
    output_file = os.path.join(config["paths"]["processed_data"], "hiv_aligned.fasta")
    
    if not os.path.exists(input_file):
        logging.error(f"Input file {input_file} not found.")
        return
        
    run_mafft(input_file, output_file, config)

if __name__ == "__main__":
    main()
