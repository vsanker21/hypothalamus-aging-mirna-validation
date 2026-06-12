"""Probe BIL RSTE3 download directory structure."""
import re
import requests

BASE = "https://download.brainimagelibrary.org/group/20240509/"
r = requests.get(BASE, timeout=60)
print("status", r.status_code, "len", len(r.text))
# parse apache directory listing
links = re.findall(r'href="([^"]+)"', r.text)
for l in sorted(set(links)):
    if l.startswith("?") or l.startswith("/"):
        continue
    print(l)
