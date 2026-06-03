"""
HIV-1 Virology Analysis: Subtyping and Drug Resistance Detection.

This module integrates with the COMET API for subtyping and the Stanford
HIVdb GraphQL API for identifying drug resistance mutations. It also
provides local fallback subtyping and generates a comprehensive HTML report.

Biological Context:
    Subtyping is essential for understanding the molecular epidemiology of
    HIV-1. Resistance mutations in the pol gene (Protease, Reverse
    Transcriptase, Integrase) determine the efficacy of antiretroviral
    therapy (ART). Monitoring these mutations is crucial for clinical
    management and public health.

Pipeline Stage:
    Phase 6 of 7: Virological Analysis.

Example:
    >>> # Run from terminal:
    >>> # python scripts/06_analysis.py
"""

from __future__ import annotations
import os
import yaml
import logging
import asyncio
import httpx
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import base64
import time
import json
import hashlib
from io import BytesIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Align
from jinja2 import Template
from typing import Any, Dict, List, Optional, Union

# Standard logging configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("logs/virology_analysis.log"),
        logging.StreamHandler()
    ]
)

class NetworkError(Exception):
    """Custom exception for network-related failures."""
    pass

class TranslationError(Exception):
    """Custom exception for protein translation failures."""
    pass

class APIRateLimitError(Exception):
    """Custom exception for API rate limit exceeding."""
    pass

def load_config() -> Dict[str, Any]:
    """
    Load pipeline configuration from a YAML file.

    Returns:
        Dict[str, Any]: Configuration dictionary.
    """
    with open("config.yaml", "r") as f:
        return yaml.safe_load(f)

def get_cache_path(sequence: Union[str, Seq], prefix: str) -> str:
    """
    Generate a deterministic cache file path based on sequence content.

    Args:
        sequence: The DNA sequence to hash.
        prefix (str): Prefix for the cache file (e.g., 'comet', 'hivdb').

    Returns:
        str: Path to the JSON cache file.
    """
    hasher = hashlib.md5()
    hasher.update(str(sequence).encode())
    return os.path.join("data/cache", f"{prefix}_{hasher.hexdigest()}.json")

async def call_comet_api(sequences: List[SeqRecord], api_url: str) -> List[Optional[Dict[str, Any]]]:
    """
    Calls COMET HIV-1 Subtyping API with caching and retry logic.

    Args:
        sequences (List[SeqRecord]): List of sequences to subtype.
        api_url (str): URL of the COMET API.

    Returns:
        List[Optional[Dict[str, Any]]]: List of subtyping results or None on failure.

    Notes:
        Ref: Struck et al. (2014) COMET: adaptive context-based modeling
        for ultra-fast HIV-1 subtyping.
    """
    results = []
    async with httpx.AsyncClient(timeout=30.0) as client:
        for record in sequences:
            cache_path = get_cache_path(record.seq, "comet")
            if os.path.exists(cache_path):
                try:
                    with open(cache_path, "r") as f:
                        results.append(json.load(f))
                    continue
                except Exception:
                    pass

            payload = {"fasta": f">{record.id}\n{str(record.seq)}"}
            success = False
            for attempt in range(3):
                try:
                    response = await client.post(api_url, json=payload)
                    if response.status_code == 429:
                        raise APIRateLimitError("Rate limit exceeded")
                    response.raise_for_status()
                    data = response.json()
                    res = {
                        'id': record.id,
                        'subtype': data.get('subtype', 'Unknown'),
                        'confidence': float(data.get('confidence', 0.0)),
                        'method_used': 'COMET API',
                        'crf_details': data.get('crf_details', '')
                    }
                    with open(cache_path, "w") as f:
                        json.dump(res, f)
                    results.append(res)
                    success = True
                    break
                except Exception as e:
                    logging.warning(f"COMET API attempt {attempt+1} failed for {record.id}: {e}")
                    await asyncio.sleep(2 ** attempt)

            if not success:
                results.append(None)
            await asyncio.sleep(0.1) # 10 req/s rate limiting
    return results

def fallback_subtyping(sequences: List[SeqRecord], ref_file: str) -> List[Dict[str, Any]]:
    """
    Fallback subtyping using local alignment against reference sequences.

    Args:
        sequences (List[SeqRecord]): Sequences to subtype.
        ref_file (str): Path to the reference FASTA file.

    Returns:
        List[Dict[str, Any]]: Subtyping results based on local alignment scores.
    """
    if not os.path.exists(ref_file):
        logging.error("Reference file for fallback not found.")
        return [ {'id': s.id, 'subtype': 'Unknown', 'confidence': 0.0, 'method_used': 'None', 'crf_details': ''} for s in sequences ]

    refs = list(SeqIO.parse(ref_file, "fasta"))
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    
    results = []
    for seq in sequences:
        best_score = -1.0
        best_subtype = "Unknown"
        for ref in refs:
            score = float(aligner.score(seq.seq, ref.seq))
            if score > best_score:
                best_score = score
                best_subtype = ref.id.split('.')[0]
        
        results.append({
            'id': seq.id,
            'subtype': best_subtype,
            'confidence': 0.5,
            'method_used': 'Local Alignment (Fallback)',
            'crf_details': 'N/A'
        })
    return results

def identify_subtypes(fasta_file: str, results_path: str) -> pd.DataFrame:
    """
    Identify HIV-1 subtypes using API calls and local fallback.

    Args:
        fasta_file (str): Path to the filtered FASTA file.
        results_path (str): Directory to save subtyping results and plots.

    Returns:
        pd.DataFrame: DataFrame containing IDs and identified subtypes.
    """
    start_time = time.time()
    logging.info("Starting HIV-1 subtype identification...")
    config = load_config()
    api_url = config['analysis']['comet_api_url']

    sequences = list(SeqIO.parse(fasta_file, "fasta"))

    try:
        loop = asyncio.get_event_loop()
        api_results = loop.run_until_complete(call_comet_api(sequences, api_url))
    except Exception as e:
        logging.warning(f"API Subtyping failed or skipped: {e}")
        api_results = [None] * len(sequences)

    final_results = []
    to_fallback = []

    for i, res in enumerate(api_results):
        if i < len(sequences):
            if res:
                final_results.append(res)
            else:
                to_fallback.append(sequences[i])

    if to_fallback:
        ref_file = "data/cache/hiv_refs.fasta"
        fallback_res = fallback_subtyping(to_fallback, ref_file)
        final_results.extend(fallback_res)

    df = pd.DataFrame(final_results)
    csv_path = os.path.join(results_path, "subtypes.csv")
    df.to_csv(csv_path, index=False)

    plt.figure(figsize=(8, 8))
    df['subtype'].value_counts().plot.pie(autopct='%1.1f%%')
    plt.title("HIV-1 Subtype Distribution")
    plt.savefig(os.path.join(results_path, "subtype_distribution.png"))
    plt.close()

    logging.info(f"Subtyping completed in {time.time() - start_time:.2f}s")
    return df

def translate_pol(sequence: str) -> str:
    """
    Translate DNA pol sequence to protein in the best reading frame.

    Args:
        sequence (str): DNA sequence string.

    Returns:
        str: Translated amino acid sequence.

    Raises:
        TranslationError: If no valid protein could be translated.
    """
    best_protein = ""
    max_length = 0
    for frame in range(3):
        dna_seq = sequence[frame:]
        trim = len(dna_seq) % 3
        if trim > 0: dna_seq = dna_seq[:-trim]
        try:
            protein = str(Seq(dna_seq).translate())
            current_len = len(protein.split('*')[0])
            if current_len > max_length:
                max_length = current_len
                best_protein = protein
        except Exception:
            continue
    if not best_protein:
        raise TranslationError("Could not translate sequence")
    return best_protein

async def call_hivdb_api(sequence_record: SeqRecord, api_url: str) -> Dict[str, Any]:
    """
    Calls Stanford HIVdb GraphQL API (Sierra) for drug resistance.

    Args:
        sequence_record (SeqRecord): Sequence to analyze.
        api_url (str): Sierra API URL.

    Returns:
        Dict[str, Any]: JSON response from the API.
    """
    query = """
    query Resistance($sequences: [UnalignedSequenceInput]!) {
      sequenceAnalysis(sequences: $sequences) {
        alignedGeneSequences {
          gene { name }
          mutations { text }
        }
        drugResistance {
          gene { name }
          drugScores {
            drug { name, displayAbbr, drugClass { name } }
            score
            level
            text
          }
        }
      }
    }
    """
    async with httpx.AsyncClient(timeout=60.0) as client:
        variables = {"sequences": [{"header": sequence_record.id, "sequence": str(sequence_record.seq)}]}
        payload = {"query": query, "variables": variables}
        response = await client.post(api_url, json=payload)
        if response.status_code == 429:
            raise APIRateLimitError("Stanford HIVdb Rate Limit")
        response.raise_for_status()
        return response.json()

def check_resistance_mutations(fasta_file: str) -> pd.DataFrame:
    """
    Detect drug resistance mutations using Stanford HIVdb API.

    Args:
        fasta_file (str): Path to the input FASTA file.

    Returns:
        pd.DataFrame: Report containing drugs, scores, and resistance levels.
    """
    start_time = time.time()
    logging.info("Starting drug resistance mutation detection...")
    config = load_config()
    api_url = config['analysis']['hivdb_api_url']
    results_path = config['paths']['results']

    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    resistance_data = []

    for record in sequences:
        cache_path = get_cache_path(record.seq, "hivdb")
        data = None
        if os.path.exists(cache_path):
            try:
                with open(cache_path, "r") as f:
                    data = json.load(f)
            except Exception:
                pass

        if data is None:
            try:
                translate_pol(str(record.seq))
                loop = asyncio.get_event_loop()
                data = loop.run_until_complete(call_hivdb_api(record, api_url))
                with open(cache_path, "w") as f:
                    json.dump(data, f)
            except Exception as e:
                logging.error(f"HIVdb analysis failed for {record.id}: {e}")
                continue

        if data and 'data' in data and data['data']['sequenceAnalysis']:
            analysis_result = data['data']['sequenceAnalysis'][0]
            # Map gene mutations
            gene_muts = {m['gene']['name']: ", ".join([mut['text'] for mut in m['mutations']])
                         for m in analysis_result['alignedGeneSequences']}

            for analysis in analysis_result['drugResistance']:
                gene = analysis['gene']['name']
                muts = gene_muts.get(gene, "")
                for score in analysis['drugScores']:
                    resistance_data.append({
                        'seq_id': record.id,
                        'gene': gene,
                        'mutation': muts,
                        'drug_class': score['drug']['drugClass']['name'],
                        'drug_name': score['drug']['name'],
                        'resistance_level': score['level'],
                        'score': score['score']
                    })

    df = pd.DataFrame(resistance_data) if resistance_data else pd.DataFrame(columns=['seq_id', 'gene', 'mutation', 'drug_class', 'drug_name', 'resistance_level', 'score'])
    csv_path = os.path.join(results_path, "resistance_report.csv")
    df.to_csv(csv_path, index=False)

    if not df.empty:
        pivot_df = df.pivot_table(index='seq_id', columns='drug_name', values='score', fill_value=0)
        plt.figure(figsize=(12, 8))
        sns.heatmap(pivot_df, annot=True, cmap="YlOrRd")
        plt.title("Stanford HIVdb Drug Resistance Scores")
        plt.tight_layout()
        plt.savefig(os.path.join(results_path, "resistance_heatmap.png"))
        plt.close()

    logging.info(f"Resistance analysis completed in {time.time() - start_time:.2f}s")
    return df

def generate_report(
    subtypes_df: pd.DataFrame,
    resistance_df: pd.DataFrame,
    results_path: str
) -> None:
    """
    Generate a consolidated HTML report with figures and tables.

    Args:
        subtypes_df (pd.DataFrame): DataFrame of subtyping results.
        resistance_df (pd.DataFrame): DataFrame of drug resistance results.
        results_path (str): Directory where the HTML report will be saved.
    """
    start_time = time.time()
    logging.info("Generating HTML summary report...")

    def img_to_base64(path: str) -> str:
        if os.path.exists(path):
            with open(path, "rb") as f:
                return base64.b64encode(f.read()).decode()
        return ""

    template_str = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>HIV-1 Analysis Report</title>
        <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.5.2/css/bootstrap.min.css">
        <style>
            body { background-color: #f8f9fa; padding: 20px; }
            .container { background-color: white; padding: 30px; border-radius: 10px; box-shadow: 0 0 10px rgba(0,0,0,0.1); }
            .img-container { text-align: center; margin: 20px 0; }
            .img-container img { max-width: 100%; height: auto; }
        </style>
    </head>
    <body>
        <div class="container">
            <h1>HIV-1 Analysis Summary</h1>
            <p>Sequences: {{ total_sequences }} | Subtypes: {{ unique_subtypes }}</p>
            <div class="img-container">
                <img src="data:image/png;base64,{{ subtype_img }}" alt="Subtypes">
            </div>
            {{ subtypes_table | safe }}
            <div class="img-container">
                <img src="data:image/png;base64,{{ resistance_img }}" alt="Resistance">
            </div>
            {{ resistance_table | safe }}
        </div>
    </body>
    </html>
    """

    template = Template(template_str)
    html_out = template.render(
        total_sequences=len(subtypes_df),
        unique_subtypes=subtypes_df['subtype'].nunique(),
        subtype_img=img_to_base64(os.path.join(results_path, "subtype_distribution.png")),
        resistance_img=img_to_base64(os.path.join(results_path, "resistance_heatmap.png")),
        subtypes_table=subtypes_df.to_html(classes="table table-sm", index=False),
        resistance_table=resistance_df.to_html(classes="table table-sm", index=False)
    )

    with open(os.path.join(results_path, "summary_report.html"), "w") as f:
        f.write(html_out)
    logging.info(f"Report generated in {time.time() - start_time:.2f}s")

def main() -> None:
    """
    Entry point for the virology analysis script.
    """
    config = load_config()
    input_file = os.path.join(config["paths"]["processed_data"], "hiv_filtered.fasta")
    results_path = config["paths"]["results"]
    if not os.path.exists(results_path):
        os.makedirs(results_path)
    subtypes_df = identify_subtypes(input_file, results_path)
    resistance_df = check_resistance_mutations(input_file)
    generate_report(subtypes_df, resistance_df, results_path)

if __name__ == "__main__":
    main()
