"""
Phylogenetic Reconstruction for HIV-1 using IQ-TREE 2.

This module performs Maximum Likelihood (ML) phylogenetic tree reconstruction
from aligned HIV-1 sequences. It includes automated model selection and
UltraFast Bootstrap (UFBoot) for node support evaluation.

Biological Context:
    Phylogeny is the cornerstone of evolutionary analysis. For HIV-1, it
    allows the identification of transmission clusters, the study of
    viral spread (phylodynamics), and the classification of subtypes.
    Maximum Likelihood methods are the gold standard for robust and
    computationally efficient tree inference.

Pipeline Stage:
    Phase 5 of 7: Phylogenetic Inference.

Example:
    >>> # Run from terminal:
    >>> # python scripts/05_phylogeny.py
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
        logging.FileHandler("logs/phylogeny.log"),
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

def run_iqtree(input_file: str, config: Dict[str, Any], output_prefix: str) -> None:
    """
    Run IQ-TREE 2 for model selection and tree reconstruction.

    Args:
        input_file (str): Path to the aligned FASTA file.
        config (Dict[str, Any]): Configuration dictionary with phylogeny parameters.
        output_prefix (str): Prefix for the generated output files.

    Notes:
        Uses `-m MFP` to select the best model and proceed with tree building.
        Employs UltraFast Bootstrap (`-B`) for evaluating branch support.
    """
    iqtree_exe = config["phylogeny"]["executable"]
    model = config["phylogeny"]["model"]
    bootstrap = config["phylogeny"]["bootstrap"]
    
    logging.info(f"Starting IQ-TREE analysis on {input_file}...")
    
    # -m MFP: ModelFinder Plus (Selects best model and proceeds)
    # -B 1000: UltraFast Bootstrap
    # -nt AUTO: Automatic thread selection
    cmd = [
        iqtree_exe,
        "-s", input_file,
        "-m", model,
        "-B", str(bootstrap),
        "-nt", "AUTO",
        "--prefix", output_prefix,
        "-redo" # Overwrite existing results
    ]
    
    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            logging.error(f"IQ-TREE error: {result.stderr}")
        else:
            logging.info(f"IQ-TREE analysis completed. Best model found and tree built.")
            logging.info(f"Results saved with prefix: {output_prefix}")
            
    except Exception as e:
        logging.error(f"Failed to run IQ-TREE: {str(e)}")

def main() -> None:
    """
    Entry point for the phylogeny script.
    """
    config = load_config()
    input_file = os.path.join(config["paths"]["processed_data"], "hiv_aligned.fasta")
    output_prefix = os.path.join(config["paths"]["results"], "hiv_phylogeny")
    
    if not os.path.exists(input_file):
        logging.error(f"Input file {input_file} not found.")
        return
        
    run_iqtree(input_file, config, output_prefix)

if __name__ == "__main__":
    main()
