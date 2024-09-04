# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

from binette import __version__

project = 'Binette'
copyright = '2024, Jean Mainguy'
author = 'Jean Mainguy'
release = __version__


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    #  
    # "sphinxcontrib.jquery",
    "sphinx.ext.duration",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.autodoc",
    'sphinx_search.extension',
     'sphinx_togglebutton',
    # "myst_nb",
    "myst_parser",
    'nbsphinx',
    'nbsphinx_link', 
    # 'sphinx.ext.napoleon',
    # 'sphinx.ext.viewcode',
    'sphinxcontrib.mermaid'
]
myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    "replacements",
    "smartquotes",
    "strikethrough",
    "substitution",
    "tasklist",
]

source_suffix = {
    '.md': 'markdown',
}


templates_path = ['_templates']

nb_execution_mode = "off"
nbsphinx_execute = 'never'
# Prefix document path to section labels, to use:
# `path/to/file:heading` instead of just `heading`
autosectionlabel_prefix_document = True

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'build', "jupyter_execute"]



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'sphinx_rtd_theme' #'alabaster' # 
html_theme = 'sphinx_rtd_theme' #'sphinx_book_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']





# Include the Plotly JavaScript in the HTML output
nbsphinx_requirejs_path = ""

# Ensures that the `require.js` is loaded for Plotly to function correctly
nbsphinx_requirejs_options = {
    'paths': {
        'plotly': 'https://cdn.plot.ly/plotly-latest.min'
    },
    'shim': {
        'plotly': {
            'exports': 'Plotly'
        }
    }
}

# Specify the default language for syntax highlighting in Sphinx
highlight_language = 'python'

# -- Options for HTML output -------------------------------------------------


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Add plotly renderer options
nbsphinx_prolog = r"""
.. raw:: html

    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
"""



