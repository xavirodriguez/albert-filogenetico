import os
import yaml
import logging
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("logs/qc.log"),
        logging.StreamHandler()
    ]
)

def load_config():
    with open("config.yaml", "r") as f:
        return yaml.safe_load(f)

def run_qc(input_file, min_len, max_ns, processed_path):
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

def generate_plots(df, output_path):
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

def main():
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
