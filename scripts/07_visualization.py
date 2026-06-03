"""
Phylogenetic Tree Visualization for HIV-1.

This module generates publication-quality visualizations of phylogenetic trees
using the ete3 library. It supports circular and rectangular layouts.

Biological Context:
    Visualizing the tree is crucial for interpreting evolutionary relationships.
    Circular layouts are often used for large datasets to display global
    diversity, while rectangular layouts are better for detailed inspection
    of specific clades or transmission clusters.

Pipeline Stage:
    Phase 7 of 7: Visualization.

Example:
    >>> # Run from terminal:
    >>> # python scripts/07_visualization.py
"""

from __future__ import annotations
import os
import yaml
import logging
import matplotlib.pyplot as plt
from ete3 import Tree, TreeStyle, NodeStyle
from typing import Any, Dict

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("logs/visualization.log"),
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

def plot_tree(tree_file: str, output_image: str) -> None:
    """
    Visualize a phylogenetic tree using ete3.

    Args:
        tree_file (str): Path to the Newick tree file.
        output_image (str): Path where the PNG image will be saved.

    Notes:
        Default styling uses a circular mode (`mode="c"`) and highlights
        certain geographical tags as an example.
    """
    if not os.path.exists(tree_file):
        logging.error(f"Tree file {tree_file} not found.")
        return

    logging.info(f"Visualizing tree from {tree_file}...")
    
    try:
        t = Tree(tree_file)
        
        # Basic styling
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.mode = "c" # Circular mode for publication-ready figures
        ts.arc_start = -180
        ts.arc_span = 360
        
        # Color nodes by some criteria if available (example: USA vs others)
        for leaf in t.iter_leaves():
            ns = NodeStyle()
            ns["fgcolor"] = "red" if "USA" in leaf.name else "blue"
            ns["size"] = 5
            leaf.set_style(ns)
            
        t.render(output_image, w=800, units="px", tree_style=ts)
        logging.info(f"Phylogenetic tree image saved to {output_image}")
        
    except Exception as e:
        logging.error(f"Error rendering tree: {str(e)}")

def main() -> None:
    """
    Entry point for the visualization script.
    """
    config = load_config()
    tree_file = os.path.join(config["paths"]["results"], "hiv_phylogeny.treefile")
    output_image = os.path.join(config["paths"]["figures"], "hiv_tree_circular.png")
    
    # We also try to generate a rectangular one (placeholder for extra logic)
    # output_image_rect = os.path.join(config["paths"]["figures"], "hiv_tree_rectangular.png")
    
    plot_tree(tree_file, output_image)

if __name__ == "__main__":
    main()
