"""
Evolutionary Model Selection for HIV-1 Phylogeny.

This module uses IQ-TREE's ModelFinder to identify the best-fit substitution
model for the aligned HIV-1 sequences before full tree reconstruction.

Biological Context:
    Selecting an appropriate evolutionary model (e.g., GTR+I+G) is vital for
    accurate phylogenetic inference. Models account for different rates of
    nucleotide substitution and rate heterogeneity across sites, which are
    pronounced in HIV-1.

Pipeline Stage:
    Phase 4 of 7: Model Selection.

Example:
    >>> # Run from terminal:
    >>> # python scripts/04_model_selection.py
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
        logging.FileHandler("logs/model_selection.log"),
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

def run_model_selection(input_file: str, config: Dict[str, Any]) -> None:
    """
    Execute ModelFinder via IQ-TREE on an aligned FASTA file.

    Args:
        input_file (str): Path to the aligned FASTA file.
        config (Dict[str, Any]): Configuration dictionary containing IQ-TREE path.

    Notes:
        This script performs only the model selection phase (`-m MF`).
        Phase 5 (Phylogeny) typically uses `MFP` to combine selection and building.
    """
    iqtree_exe = config["phylogeny"]["executable"]
    
    logging.info(f"Running ModelTest-NG/ModelFinder via IQ-TREE on {input_file}...")
    
    cmd = [
        iqtree_exe,
        "-s", input_file,
        "-m", "MF", # ModelFinder only
        "-nt", "AUTO"
    ]
    
    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            logging.error(f"Model selection error: {result.stderr}")
        else:
            logging.info("Model selection completed.")
            # Search for best model in stdout or .iqtree file
            for line in result.stdout.split('\n'):
                if "Best-fit model:" in line:
                    logging.info(f"FOUND: {line}")
                    
    except Exception as e:
        logging.error(f"Failed to run model selection: {str(e)}")

def main() -> None:
    """
    Entry point for the model selection script.
    """
    config = load_config()
    input_file = os.path.join(config["paths"]["processed_data"], "hiv_aligned.fasta")
    
    if not os.path.exists(input_file):
        logging.error(f"Input file {input_file} not found.")
        return
        
    run_model_selection(input_file, config)

if __name__ == "__main__":
    main()
