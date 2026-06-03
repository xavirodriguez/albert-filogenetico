# IQ-TREE handles model selection automatically with -m MFP. 
# This script is a wrapper for Phase 4 specifically if only model selection is needed.

import os
import yaml
import subprocess
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("logs/model_selection.log"),
        logging.StreamHandler()
    ]
)

def load_config():
    with open("config.yaml", "r") as f:
        return yaml.safe_load(f)

def run_model_selection(input_file, config):
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

def main():
    config = load_config()
    input_file = os.path.join(config["paths"]["processed_data"], "hiv_aligned.fasta")
    
    if not os.path.exists(input_file):
        logging.error(f"Input file {input_file} not found.")
        return
        
    run_model_selection(input_file, config)

if __name__ == "__main__":
    main()
