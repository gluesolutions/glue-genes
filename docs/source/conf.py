# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "glue genes"
copyright = "2023, glue solutions"
# release = "1.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx_automodapi.automodapi",
    "sphinx.ext.intersphinx",
    "sphinx_design",
    "sphinx_rtd_theme",
    "numpydoc",
]
numpydoc_show_class_members = False

templates_path = ["_templates"]
exclude_patterns = []

# The master toctree document.
master_doc = "index"


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_logo = "images/glue_genes_icon.png"
html_theme = "sphinx_rtd_theme"
html_favicon = "images/favicon.png"

html_context = {
    "github_user": "gluesolutions",
    "github_repo": "glue-genes",
    "github_version": "main",
    "doc_path": "docs",
}

html_static_path = ["_static"]
html_css_files = ["style.css"]


# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
# html_theme_options = dict(
#    repository_url="https://github.com/gluesolutions/glue-genes",
#    repository_branch="main",
#    path_to_docs="docs",
#    use_edit_page_button=True,
#    use_repository_button=True,
#    use_issues_button=True,
#    home_page_in_toc=False,
#    extra_navbar="",
#    navbar_footer_text="",
#    extra_footer="""<br>
#    Theme by the <a href="https://ebp.jupyterbook.org">Executable Book Project</a></p>""",
# )
