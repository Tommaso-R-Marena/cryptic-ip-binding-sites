"""Sphinx configuration for cryptic_ip API docs."""

from datetime import datetime
import os
import sys

sys.path.insert(0, os.path.abspath('../..'))

project = 'cryptic-ip-binding-sites'
author = 'Tommaso R. Marena'
copyright = f"{datetime.now():%Y}, {author}"
release = '0.1.0'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
]

autosummary_generate = True
autodoc_default_options = {
    'members': True,
    'undoc-members': False,
    'show-inheritance': True,
}

html_theme = 'alabaster'
master_doc = 'index'
