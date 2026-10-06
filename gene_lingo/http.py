"""Minimal JSON HTTP helper (stdlib only). Tests replace `request_json`."""
import json
import time
import urllib.error
import urllib.parse
import urllib.request

USER_AGENT = "gene-lingo/0.1"


def request_json(url, params=None, data=None, retries=3):
    """GET (or form-POST when `data` is given) and return parsed JSON."""
    if params:
        url += ("&" if "?" in url else "?") + urllib.parse.urlencode(params)
    body = urllib.parse.urlencode(data).encode() if data is not None else None
    req = urllib.request.Request(
        url, data=body, headers={"Accept": "application/json", "User-Agent": USER_AGENT}
    )
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.load(resp)
        except urllib.error.HTTPError as e:
            if e.code == 429 and attempt < retries - 1:  # rate limited
                time.sleep(float(e.headers.get("Retry-After", 2 ** attempt)))
                continue
            if e.code == 400:  # Ensembl: unknown symbol / no homologs
                return None
            raise
