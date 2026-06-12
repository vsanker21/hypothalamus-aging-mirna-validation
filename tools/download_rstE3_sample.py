"""Download RSTE3 metadata header sample from BIL."""
import requests
from pathlib import Path

url = "https://download.brainimagelibrary.org/group/20240509/RSTE3/RSTE3_metadata.csv"
dest = Path("data/references/bil_rstE3/RSTE3_metadata.csv")
dest.parent.mkdir(parents=True, exist_ok=True)
if dest.exists() and dest.stat().st_size > 1_000_000:
    print("exists", dest, dest.stat().st_size)
else:
    print("downloading metadata...")
    with requests.get(url, stream=True, timeout=600) as r:
        r.raise_for_status()
        with open(dest, "wb") as f:
            for chunk in r.iter_content(1 << 20):
                if chunk:
                    f.write(chunk)
    print("done", dest.stat().st_size)

import pandas as pd
df = pd.read_csv(dest, nrows=5)
print("cols", list(df.columns))
print(df.head(2).T)
