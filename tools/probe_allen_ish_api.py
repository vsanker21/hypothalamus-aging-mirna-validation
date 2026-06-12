"""Probe Allen ISH RMA API — find working StructureUnionize query."""
import requests

RMA = "http://api.brain-map.org/api/v2/data/query.json"
CSV = "http://api.brain-map.org/api/v2/data/query.csv"


def q_json(criteria: str):
    r = requests.get(RMA, params={"criteria": criteria, "num_rows": "all"}, timeout=90)
    r.raise_for_status()
    msg = r.json().get("msg", [])
    if isinstance(msg, str):
        print("ERR", msg[:200])
        return []
    print("n=", len(msg))
    if msg:
        print("keys", list(msg[0].keys())[:12])
        print("sample", {k: msg[0].get(k) for k in ("id", "structure_id", "expression_energy", "expression_density") if k in msg[0]})
    return msg


def q_csv(criteria: str, tabular: str = ""):
    params = {"criteria": criteria, "num_rows": "all"}
    if tabular:
        params["tabular"] = tabular
    r = requests.get(CSV, params=params, timeout=90)
    print("csv status", r.status_code, "len", len(r.text))
    print(r.text[:800])
    return r.text


# Gpr50 mouse gene id
genes = q_json("model::Gene,rma::criteria,[acronym$eq'Gpr50'],[organism_id$eq2]")
gid = genes[0]["id"]
print("Gpr50 id", gid)

# SectionDataSet via genes association (Mouse Brain product id=1)
print("\n--- SectionDataSet via genes ---")
sds = q_json(f"model::SectionDataSet,rma::criteria,genes[id$eq{gid}],products[id$eq1]")
if sds:
    sds_id = sds[0]["id"]
    print("first SDS", sds_id, "plane", sds[0].get("plane_of_section_id"))
    print("\n--- StructureUnionize for SDS + structure 135 ---")
    q_json(
        f"model::StructureUnionize,rma::criteria,[section_data_set_id$eq{sds_id}],[structure_id$eq135]"
    )
    print("\n--- StructureUnionize via section_data_set include ---")
    q_json(
        f"model::StructureUnionize,rma::criteria,section_data_set[id$eq{sds_id}],[structure_id$eq135]"
    )

print("\n--- CSV tabular all structures for gene ---")
q_csv(
    f"model::StructureUnionize,rma::criteria,section_data_set(genes[id$eq{gid}],products[id$eq1])",
    tabular="structures.acronym,structures.name,structure_unionizes.expression_energy,structure_unionizes.structure_id",
)

print("\n--- Foxj1 ---")
fox = q_json("model::Gene,rma::criteria,[acronym$eq'Foxj1'],[organism_id$eq2]")
if fox:
    fgid = fox[0]["id"]
    q_csv(
        f"model::StructureUnionize,rma::criteria,section_data_set(genes[id$eq{fgid}],products[id$eq1])",
        tabular="structures.acronym,structures.name,structure_unionizes.expression_energy",
    )
