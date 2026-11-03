"""The app's page registry: file, sidebar title, icon and section of every page.

app.py builds the navigation from SECTIONS; pages use link() to point at each other, so titles live in one place.
"""
import streamlit as st

PAGES = {
    "home": ("views/home.py", "Start here", ":material/home:"),
    "species_explained": ("views/species_explained.py", "How species matching works", ":material/school:"),
    "species_lookup": ("views/species_lookup.py", "Look up genes across species", ":material/search:"),
    "species_custom": ("views/species_custom.py", "Compare your own genomes", ":material/biotech:"),
    "species_osdr": ("views/species_osdr.py", "Match an OSDR study", ":material/science:"),
    "strains_explained": ("views/strains_explained.py", "How strain matching works", ":material/school:"),
    "strains_compare": ("views/strains_compare.py", "Compare two strains", ":material/compare_arrows:"),
    "converter": ("views/converter.py", "Convert gene IDs", ":material/swap_horiz:"),
    "converter_how": ("views/converter_how_it_works.py", "How the converter works", ":material/school:"),
}
SECTIONS = {
    "Start": ["home"],
    "Bacteria: between species": ["species_explained", "species_lookup", "species_custom", "species_osdr"],
    "Bacteria: between strains": ["strains_explained", "strains_compare"],
    "Human & animal genes": ["converter", "converter_how"],
}


def title(key):
    return PAGES[key][1]


def pages():
    """{section: [st.Page]} for st.navigation."""
    return {sec: [st.Page(PAGES[k][0], title=PAGES[k][1], icon=PAGES[k][2], default=(k == "home")) for k in keys]
            for sec, keys in SECTIONS.items()}


def link(container, key, label=None):
    """Page link that degrades to plain text when a page runs on its own (tests) instead of through app.py."""
    path, name, icon = PAGES[key]
    try:
        container.page_link(path, label=label or name, icon=icon)
    except st.errors.StreamlitPageNotFoundError:
        container.write(label or name)
