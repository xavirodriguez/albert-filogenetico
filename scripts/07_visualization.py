import os
import yaml
import logging
import matplotlib.pyplot as plt
from ete3 import Tree, TreeStyle, NodeStyle

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("logs/visualization.log"),
        logging.StreamHandler()
    ]
)

def load_config():
    with open("config.yaml", "r") as f:
        return yaml.safe_load(f)

def plot_tree(tree_file, output_image):
    """
    Visualizes the phylogenetic tree using ete3.
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
        
        # Color nodes by some criteria if available (placeholder)
        for leaf in t.iter_leaves():
            ns = NodeStyle()
            ns["fgcolor"] = "red" if "USA" in leaf.name else "blue"
            ns["size"] = 5
            leaf.set_style(ns)
            
        t.render(output_image, w=800, units="px", tree_style=ts)
        logging.info(f"Phylogenetic tree image saved to {output_image}")
        
    except Exception as e:
        logging.error(f"Error rendering tree: {str(e)}")

def main():
    config = load_config()
    tree_file = os.path.join(config["paths"]["results"], "hiv_phylogeny.treefile")
    output_image = os.path.join(config["paths"]["figures"], "hiv_tree_circular.png")
    
    # We also try to generate a rectangular one
    output_image_rect = os.path.join(config["paths"]["figures"], "hiv_tree_rectangular.png")
    
    plot_tree(tree_file, output_image)
    
    # Additional logic for other visualizations (e.g. mutation heatmaps) could go here

if __name__ == "__main__":
    main()
