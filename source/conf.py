# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "Brightway training CEA"
copyright = "2025, V. Maneval, M. Boutrouelle"
author = "V. Maneval, M. Boutrouelle"
release = "0.1"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
]

templates_path = ["_templates"]
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "furo"
html_title = "Brightway training"
html_static_path = ["https://vmvl-cea.github.io/CEA_brightway_training/_static"]
html_theme_options = {
    "light_logo": "BW_black_transparent_square.png",
    "dark_logo": "BW_white_transparent_square.png",
}
html_context = {
    "license_name": "Creative Commons Attribution 4.0 International",
    "license_url": "https://creativecommons.org/licenses/by/4.0/",
}

html_baseurl = "https://vmvl-cea.github.io/CEA_brightway_training/"
