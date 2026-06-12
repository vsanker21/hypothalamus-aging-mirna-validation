import requests

RMA = "http://api.brain-map.org/api/v2/data/query.json"
MOUSE = 2
PROD = 1

for ac in ["Rax", "Tm4sf1", "Aqp4", "Gpr50"]:
    rows = requests.get(
        RMA,
        params={"criteria": f"model::Gene,rma::criteria,[acronym$eq'{ac}'],[organism_id$eq{MOUSE}]", "num_rows": "all"},
        timeout=60,
    ).json()["msg"]
    print(ac, "n genes", len(rows))
    for r in rows:
        gid = r["id"]
        sds = requests.get(
            RMA,
            params={"criteria": f"model::SectionDataSet,rma::criteria,genes[id$eq{gid}],products[id$eq{PROD}]", "num_rows": "all"},
            timeout=60,
        ).json()["msg"]
        n = len(sds) if isinstance(sds, list) else 0
        print(" ", gid, r.get("name"), "sds", n)
