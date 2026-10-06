from gene_lingo import convert, orthologs, http


def fake(url, params=None, data=None, retries=3):
    if "mygene.info" in url:
        q = data["q"].split(",")
        if data["species"] == 9606:
            db = {"TP53": {"symbol": "TP53", "entrezgene": "7157", "ensembl": {"gene": "ENSG00000141510"}}}
            return [dict(db[g], query=g) if g in db else {"query": g, "notfound": True} for g in q]
        return [{"query": g, "symbol": "Trp53"} for g in q]
    if "homology/id/homo_sapiens/ENSG00000141510" in url:
        assert params["target_species"] == "mus_musculus"
        return {"data": [{"homologies": [
            {"type": "ortholog_one2one",
             "target": {"id": "ENSMUSG00000059552", "species": "mus_musculus", "perc_id": 80, "perc_pos": 85}}]}]}
    raise AssertionError(url)


def test_convert(monkeypatch):
    monkeypatch.setattr(http, "request_json", fake)
    assert convert(["TP53", "NOPE"], to="entrez") == {"TP53": ["7157"], "NOPE": []}
    assert convert(["TP53"], to=["symbol", "ensembl"])["TP53"]["ensembl"] == ["ENSG00000141510"]


def test_orthologs(monkeypatch):
    monkeypatch.setattr(http, "request_json", fake)
    res = orthologs(["TP53"], "human", "mouse")
    assert res["TP53"][0]["id"] == "ENSMUSG00000059552"
    assert res["TP53"][0]["symbol"] == "Trp53"
    assert orthologs(["TP53"], "human", "mouse", "one2many")["TP53"] == []


def test_translate(monkeypatch):
    from gene_lingo.pipeline import translate
    monkeypatch.setattr(http, "request_json", fake)
    r = translate(["TP53"], "human", "mouse", "symbol", "symbol")
    assert r[0]["output_id"] == "Trp53" and r[0]["output_species"] == "mus_musculus"
    r = translate(["TP53", "NOPE"], "human", "human", "auto", "entrez")
    assert [x["output_id"] for x in r] == ["7157", ""]
