07_visualization
================

.. py:module:: 07_visualization

.. autoapi-nested-parse::

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

   .. admonition:: Example

      >>> # Run from terminal:
      >>> # python scripts/07_visualization.py



Functions
---------

.. autoapisummary::

   07_visualization.load_config
   07_visualization.plot_tree
   07_visualization.main


Module Contents
---------------

.. py:function:: load_config() -> Dict[str, Any]

   Load pipeline configuration from a YAML file.

   :returns: Configuration dictionary.
   :rtype: Dict[str, Any]


.. py:function:: plot_tree(tree_file: str, output_image: str) -> None

   Visualize a phylogenetic tree using ete3.

   :param tree_file: Path to the Newick tree file.
   :type tree_file: :py:class:`str`
   :param output_image: Path where the PNG image will be saved.
   :type output_image: :py:class:`str`

   .. admonition:: Notes

      Default styling uses a circular mode (`mode="c"`) and highlights
      certain geographical tags as an example.


.. py:function:: main() -> None

   Entry point for the visualization script.
