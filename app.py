"""Gene Lingo: Streamlit entry point.  streamlit run app.py

Only the navigation lives here; the pages are in views/ and are registered in gene_lingo/nav.py.
"""
import streamlit as st

from gene_lingo import nav

st.set_page_config(page_title="Gene Lingo", layout="wide")
st.navigation(nav.pages()).run()
