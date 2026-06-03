import os
import yaml
import subprocess
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("logs/phylogeny.log"),
        logging.StreamHandler()
    ]
)

def load_config():
    with open("config.yaml", "r") as f:
        return yaml.safe_load(f)

def run_iqtree(input_file, config, output_prefix):
    """
    Runs IQ-TREE 2 for model selection and tree reconstruction.
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

def main():
    config = load_config()
    input_file = os.path.join(config["paths"]["processed_data"], "hiv_aligned.fasta")
    output_prefix = os.path.join(config["paths"]["results"], "hiv_phylogeny")
    
    if not os.path.exists(input_file):
        logging.error(f"Input file {input_file} not found.")
        return
        
    run_iqtree(input_file, config, output_prefix)

if __name__ == "__main__":
    main()
