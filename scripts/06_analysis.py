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

# Standard logging configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("logs/virology_analysis.log"),
        logging.StreamHandler()
    ]
)

class NetworkError(Exception): pass
class TranslationError(Exception): pass
class APIRateLimitError(Exception): pass

def load_config():
    with open("config.yaml", "r") as f:
        return yaml.safe_load(f)

def get_cache_path(sequence, prefix):
    hasher = hashlib.md5()
    hasher.update(str(sequence).encode())
    return os.path.join("data/cache", f"{prefix}_{hasher.hexdigest()}.json")

async def call_comet_api(sequences, api_url):
    """
    Calls COMET HIV-1 Subtyping API.
    Ref: Struck et al. (2014) COMET: adaptive context-based modeling for ultra-fast HIV-1 subtyping.
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
                except: pass

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
            await asyncio.sleep(0.1) # 10 req/s
    return results

def fallback_subtyping(sequences, ref_file):
    """
    Fallback subtyping using local alignment against reference sequences.
    """
    if not os.path.exists(ref_file):
        logging.error("Reference file for fallback not found.")
        return [ {'id': s.id, 'subtype': 'Unknown', 'confidence': 0.0, 'method_used': 'None', 'crf_details': ''} for s in sequences ]

    refs = list(SeqIO.parse(ref_file, "fasta"))
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    
    results = []
    for seq in sequences:
        best_score = -1
        best_subtype = "Unknown"
        for ref in refs:
            score = aligner.score(seq.seq, ref.seq)
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

def identify_subtypes(fasta_file, results_path):
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

def translate_pol(sequence):
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
        except Exception: continue
    if not best_protein: raise TranslationError("Could not translate sequence")
    return best_protein

async def call_hivdb_api(sequence_record, api_url):
    """
    Calls Stanford HIVdb GraphQL API (Sierra).
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

def check_resistance_mutations(fasta_file):
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
            except: pass

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

def generate_report(subtypes_df, resistance_df, results_path):
    start_time = time.time()
    logging.info("Generating HTML summary report...")

    def img_to_base64(path):
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

def main():
    config = load_config()
    input_file = os.path.join(config["paths"]["processed_data"], "hiv_filtered.fasta")
    results_path = config["paths"]["results"]
    if not os.path.exists(results_path): os.makedirs(results_path)
    subtypes_df = identify_subtypes(input_file, results_path)
    resistance_df = check_resistance_mutations(input_file)
    generate_report(subtypes_df, resistance_df, results_path)

if __name__ == "__main__":
    main()
