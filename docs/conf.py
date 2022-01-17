# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
master_doc = 'index'
# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'nbis-meta'
copyright = '2022, John Sundh'
author = 'John Sundh'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
#extensions = [
#    'sphinx_rtd_theme'
#]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_logo = 'https://nbis.se/assets/img/logos/nbislogo-green-txt.svg'

html_theme_options = {
    #'rightsidebar': False,
    #'stickysidebar': True,
    #'collapsiblesidebar': True,
    #'textcolor': '#000000',
    #'linkcolor': '#f47d20',
    # sidebar
    #'sidebarbgcolor': '#EBEBEB',
    #'sidebartextcolor': '#000000',
    #'sidebarbtncolor': '#D6D6D6',
    #'sidebarlinkcolor': '#f47d20',
    'sidebarwidth': '25%',
    # footer
    #'footerbgcolor': '#EBEBEB',
    #'footertextcolor': '#000000',
    # relbar
    #'relbarbgcolor': '#EBEBEB',
    #'relbartextcolor': '#000000',
    #'relbarlinkcolor': '#f47d20',
    # headers
    #'headbgcolor': '#FFFFFF',

}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

root_doc = "index"

html_sidebars = {
   '**': ['globaltoc.html', 'sourcelink.html', 'searchbox.html'],
    'index': ['globaltoc.html', 'localtoc.html', 'sourcelink.html', 'searchbox.html'],
   'using/windows': ['windowssidebar.html', 'searchbox.html'],
}
