# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys
import os
sys.path.insert(0, os.path.abspath(".."))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Geoclide'
copyright = '2025, Geoclide team'
author = 'Geoclide team'
release = '2.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.todo',
              'sphinx.ext.autodoc',
              'myst_parser',
              'nbsphinx',
              "sphinx.ext.napoleon"]

templates_path = ['_templates']

# autodoc_default_flags = ['members', 'show-inheritance', 'special-members']

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '*constante.rst',
                    'modules.rst']

autodoc_default_options = {
    'members': True,
    'show-inheritance': True,
    'special-members': '__call__'
}

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}

numpydoc_show_class_members = True

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
