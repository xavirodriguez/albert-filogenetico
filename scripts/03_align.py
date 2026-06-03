import os
import yaml
import subprocess
import logging
from Bio import SeqIO

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("logs/alignment.log"),
        logging.StreamHandler()
    ]
)

def load_config():
    with open("config.yaml", "r") as f:
        return yaml.safe_load(f)

def run_mafft(input_file, output_file, config):
    """
    Executes MAFFT alignment with configured parameters.
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

def main():
    config = load_config()
    input_file = os.path.join(config["paths"]["processed_data"], "hiv_filtered.fasta")
    output_file = os.path.join(config["paths"]["processed_data"], "hiv_aligned.fasta")
    
    if not os.path.exists(input_file):
        logging.error(f"Input file {input_file} not found.")
        return
        
    run_mafft(input_file, output_file, config)

if __name__ == "__main__":
    main()
