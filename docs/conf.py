# === FILE: docs/conf.py ===
import os
import sys
from datetime import datetime

# Path setup
sys.path.insert(0, os.path.abspath('..'))

# Project information
project = "HIV-1 Genomic Pipeline"
author = "Xavier Rodríguez"
copyright = f"{datetime.now().year}, {author}"
release = "1.0.0"

# General configuration
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.doctest",
    "myst_parser",
    "autoapi.extension",
    "sphinx_copybutton",
    "sphinx_autodoc_typehints",
    "sphinxcontrib.mermaid",
]

# AutoAPI configuration
autoapi_dirs = ["../scripts"]
autoapi_options = [
    "members",
    "undoc-members",
    "private-members",
    "show-inheritance",
    "show-module-summary",
    "special-members",
    "imported-members",
]
autoapi_python_class_content = "both"
autoapi_add_toctree_entry = False
autoapi_keep_files = True

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = True
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = True
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = True

# Intersphinx mapping
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
}

# MyST settings
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "html_image",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# HTML output options
html_theme = "furo"
html_static_path = ["_static"]
html_css_files = ["custom.css"]

# Source suffixes
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# Ignore warnings that are not critical for verification
suppress_warnings = ["toc.no_title", "intersphinx.inventory"]
