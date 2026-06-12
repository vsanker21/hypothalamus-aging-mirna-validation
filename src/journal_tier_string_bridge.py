"""
STRING protein–protein association graph bridge (mouse).

Fetches a compact STRING network (public API) for:
  - Curated Piezo1 / mechanistic priors (Ca²⁺, ERK, YAP/TAZ, NF-κB, TGF-β, HIF axes; mouse symbols)
  - Top miRNA-target–union genes by weighted burden (from exploratory_crossmodal CSV)

Primary metrics (exploratory):
  (1) n_high_conf_edges — edges with STRING combined score >= score_min between the
      mechanism seed set and the target-union gene subset (both must appear in the returned network).
  (2) stress_proximity — median graph distance from a stress-pathway gene set (PROGENy-derived
      pathway names matching TNF/NFκB/TGF/HIF/apoptosis/JAK patterns) to the Piezo1 hub
      (Piezo1 ∪ its 1-hop neighbors in the induced graph), restricted to genes in the network.

Null for (1): when the node pool is large enough, uniform random samples **without replacement**
of the same size as the high-burden union subset (excluding fixed mechanism seeds); otherwise
degree-bin–matched draws with replacement (may include duplicate genes in the control set).

Requires network access to string-db.org unless STRING_BRIDGE_OFFLINE=1 (writes skip JSON).

Caveat: STRING is functional association, not causal directionality; mouse ID mapping is via STRING.
"""
from __future__ import annotations

import json
import math
import os
import re
from collections import defaultdict, deque
from pathlib import Path

import numpy as np
import pandas as pd
import requests


STRING_NETWORK_URL = "https://string-db.org/api/tsv/network"

# Literature-aligned priors (mouse symbols; STRING resolves aliases).
MECHANISM_SEEDS: tuple[str, ...] = (
    "Piezo1",
    "Yap1",
    "Wwtr1",
    "Nfkbia",
    "Rela",
    "Nfkb1",
    "Smad2",
    "Smad3",
    "Smad4",
    "Tgfbr1",
    "Tgfbr2",
    "Hif1a",
    "Epas1",
    "Trpv4",
    "Calm1",
    "Camk2d",
    "Mapk1",
    "Mapk3",
    "Ppp3ca",
    "Nfatc1",
    "Ccn1",
    "Ccn2",
)

PROG_STRESS_RE = re.compile(
    r"(?i)(tgfb|tnf|nfkb|nf-?kb|hypox|hif|il-?6|interferon|apopt|p53|jak|stat3|inflam|toll|tlr)"
)


def _fetch_string_network(identifiers: list[str], *, required_score: int = 400, timeout: int = 180) -> pd.DataFrame:
    """POST identifiers (mouse, NCBI taxon 10090)."""
    payload = {
        "identifiers": "\r".join(identifiers),
        "species": "10090",
        "network_type": "functional",
        "required_score": str(int(required_score)),
        "add_nodes": "0",
    }
    r = requests.post(STRING_NETWORK_URL, data=payload, timeout=timeout)
    r.raise_for_status()
    lines = [ln for ln in r.text.strip().splitlines() if ln.strip()]
    if len(lines) < 2:
        return pd.DataFrame()
    header = lines[0].split("\t")
    rows = []
    for ln in lines[1:]:
        parts = ln.split("\t")
        if len(parts) < len(header):
            continue
        rows.append(dict(zip(header, parts)))
    return pd.DataFrame(rows)


def _norm_sym(s: str) -> str:
    return str(s).strip().upper()


def _build_graph(df: pd.DataFrame, score_col: str = "score", score_min: float = 0.35) -> dict[str, list[tuple[str, float]]]:
    adj: dict[str, list[tuple[str, float]]] = defaultdict(list)
    if df.empty or score_col not in df.columns:
        return adj
    c_a = "preferredName_A" if "preferredName_A" in df.columns else None
    c_b = "preferredName_B" if "preferredName_B" in df.columns else None
    if not c_a or not c_b:
        return adj
    for _, row in df.iterrows():
        try:
            sc = float(row[score_col])
        except (TypeError, ValueError):
            continue
        if sc < score_min:
            continue
        a = _norm_sym(row[c_a])
        b = _norm_sym(row[c_b])
        if not a or not b or a == b:
            continue
        adj[a].append((b, sc))
        adj[b].append((a, sc))
    return adj


def _degree(adj: dict[str, list[tuple[str, float]]]) -> dict[str, int]:
    return {n: len(adj.get(n, ())) for n in adj}


def _bfs_distances(adj: dict[str, list[tuple[str, float]]], sources: set[str]) -> dict[str, int]:
    dist: dict[str, int] = {}
    q = deque()
    for s in sources:
        if s in adj:
            dist[s] = 0
            q.append(s)
    while q:
        u = q.popleft()
        du = dist[u]
        for v, _ in adj.get(u, ()):
            if v not in dist:
                dist[v] = du + 1
                q.append(v)
    return dist


def _stress_genes_from_progeny(burden_genes: set[str]) -> set[str]:
    try:
        import decoupler as dc
    except ImportError:
        return set()
    net = dc.op.progeny(organism="mouse", license="academic", verbose=False)
    out: set[str] = set()
    for pw, sub in net.groupby("source", sort=False):
        if not PROG_STRESS_RE.search(str(pw)):
            continue
        for t in sub["target"]:
            g = _norm_sym(t)
            if g in burden_genes:
                out.add(g)
    return out


def _deg_bin(d: int) -> int:
    return int(round(math.log2(max(d, 0) + 1)))


def _degree_matched_draw(
    pool: list[str],
    deg_map: dict[str, int],
    targets: list[str],
    rng: np.random.Generator,
) -> list[str]:
    """For each gene in targets, pick a random gene from pool with same degree bin (fallback uniform)."""
    bins: dict[int, list[str]] = defaultdict(list)
    for g in pool:
        bins[_deg_bin(deg_map.get(g, 0))].append(g)
    draws: list[str] = []
    for g in targets:
        b = _deg_bin(deg_map.get(g, 0))
        cand = [x for x in bins.get(b, ()) if x != g]
        if not cand:
            cand = [x for x in pool if x != g]
        if not cand:
            continue
        draws.append(rng.choice(cand))
    return draws


def _n_edges_between(adj: dict[str, list[tuple[str, float]]], set_a: set[str], set_b: set[str]) -> int:
    seen = set()
    n = 0
    for u in set_a:
        if u not in adj:
            continue
        for v, _ in adj[u]:
            if v in set_b:
                e = (u, v) if u <= v else (v, u)
                if e not in seen:
                    seen.add(e)
                    n += 1
    return n


def run_string_piezo_bridge(
    out_dir: Path,
    log,
    *,
    top_union_genes: int = 380,
    n_null: int = 400,
    score_min: float = 0.35,
    string_required_score: int = 400,
    seed: int = 46,
) -> None:
    log("\n=== STRING graph bridge (Piezo1 priors × miRNA target union) ===")
    if os.environ.get("STRING_BRIDGE_OFFLINE", "").strip().lower() in ("1", "true", "yes"):
        log("STRING bridge skipped (STRING_BRIDGE_OFFLINE=1).")
        return

    bur_path = out_dir / "exploratory_crossmodal_gene_burden_vs_aging_logfc.csv"
    uni_path = out_dir / "mirna_target_union_genes.csv"
    if not bur_path.is_file() or not uni_path.is_file():
        log("STRING bridge: missing burden CSV or mirna_target_union_genes.csv; skipped.")
        return

    bdf = pd.read_csv(bur_path)
    bdf["gene_u"] = bdf["gene"].astype(str).str.strip().str.upper()
    bdf = bdf.sort_values("weighted_burden", ascending=False)
    top_genes = bdf["gene_u"].head(top_union_genes).tolist()
    burden_rank = dict(zip(bdf["gene_u"], bdf["weighted_burden"].astype(float)))

    uni = pd.read_csv(uni_path)
    col = "gene" if "gene" in uni.columns else uni.columns[0]
    union_all = { _norm_sym(x) for x in uni[col].astype(str) }

    seeds = [_norm_sym(s) for s in MECHANISM_SEEDS]
    seeds_u = list(dict.fromkeys(seeds))

    id_list = list(dict.fromkeys(seeds_u + top_genes))
    if len(id_list) < 10:
        log("STRING bridge: too few identifiers; skipped.")
        return

    try:
        df = _fetch_string_network(id_list, required_score=string_required_score)
    except Exception as e:
        log(f"STRING bridge: network request failed ({e}); skipped.")
        return

    if df.empty:
        log("STRING bridge: empty network response; skipped.")
        return

    adj = _build_graph(df, score_min=score_min)
    nodes = set(adj.keys())
    if not nodes:
        log("STRING bridge: no edges after score filter; skipped.")
        return

    deg_map = _degree(adj)
    union_sub = set(top_genes) & nodes & union_all
    seed_in = set(seeds_u) & nodes
    pool = [n for n in nodes if n not in seed_in]
    trimmed = False
    if len(union_sub) > len(pool) and pool:
        us = sorted(union_sub, key=lambda g: -float(burden_rank.get(g, 0.0)))
        union_sub = set(us[: len(pool)])
        trimmed = True

    pz = _norm_sym("Piezo1")
    if pz not in nodes:
        log("STRING bridge: Piezo1 absent from induced network after API filter; hub distances may be empty.")

    hub: set[str] = set()
    if pz in nodes:
        hub.add(pz)
        for v, _ in adj.get(pz, ()):
            hub.add(v)

    stress_raw = _stress_genes_from_progeny(set(bdf["gene_u"]))
    stress = stress_raw & nodes

    dist_from_hub: dict[str, int] = {}
    if hub:
        dist_from_hub = _bfs_distances(adj, hub)

    obs_edges = _n_edges_between(adj, seed_in, union_sub)
    stress_dists = [dist_from_hub[g] for g in stress if g in dist_from_hub]
    obs_med_dist = float(np.median(stress_dists)) if len(stress_dists) >= 5 else float("nan")

    rng = np.random.default_rng(seed)
    union_list = sorted(union_sub)
    null_edges: list[float] = []
    if len(pool) >= len(union_list) and union_list:
        for _ in range(n_null):
            fake_arr = rng.choice(np.array(pool, dtype=object), size=len(union_list), replace=False)
            fake_set = set(fake_arr.tolist())
            null_edges.append(float(_n_edges_between(adj, seed_in, fake_set)))
        arr_ne = np.array(null_edges, dtype=float)
        perm_p_edges = float(np.mean(arr_ne >= obs_edges)) if arr_ne.size else float("nan")
        edge_mode = "uniform_without_replacement"
    elif union_list:
        for _ in range(n_null):
            fake = _degree_matched_draw(pool, deg_map, union_list, rng)
            fake_set = set(fake)
            null_edges.append(float(_n_edges_between(adj, seed_in, fake_set)))
        arr_ne = np.array(null_edges, dtype=float)
        perm_p_edges = float(np.mean(arr_ne >= obs_edges)) if arr_ne.size else float("nan")
        edge_mode = "degree_matched_with_replacement_draw"
    else:
        perm_p_edges = float("nan")
        edge_mode = "none"

    # Stress-distance null: random gene sets (same |stress|) drawn from pool \ hub vs hub distances
    null_med = []
    perm_p_dist = float("nan")
    if len(stress) >= 5 and dist_from_hub and pool:
        stress_k = len(stress)
        pool_no_hub = [g for g in pool if g not in hub]
        m_draw = min(stress_k, len(pool_no_hub))
        if m_draw >= 5:
            for _ in range(n_null):
                draw = rng.choice(pool_no_hub, size=m_draw, replace=False)
                dlist = [dist_from_hub.get(g) for g in draw if g in dist_from_hub]
                dlist = [x for x in dlist if x is not None]
                if len(dlist) >= 5:
                    null_med.append(float(np.median(dlist)))
            if null_med and np.isfinite(obs_med_dist):
                perm_p_dist = float(np.mean(np.array(null_med) <= obs_med_dist))

    summ = {
        "n_string_nodes_induced": len(nodes),
        "n_string_edges_passing_filter": sum(len(v) for v in adj.values()) // 2,
        "score_min": score_min,
        "string_required_score_param": string_required_score,
        "n_union_genes_in_network": len(union_sub),
        "union_subset_trimmed_to_pool_size": trimmed,
        "n_mechanism_seeds_in_network": len(seed_in),
        "n_edges_mechanism_seeds_to_union_subset": obs_edges,
        "perm_p_edges_ge_obs_uniform_sample_null": perm_p_edges,
        "edge_null_mode": edge_mode,
        "n_stress_progeny_genes_in_network": len(stress),
        "median_graph_dist_stress_to_piezo1_hub": obs_med_dist,
        "perm_p_median_dist_le_obs_stress_label_null": perm_p_dist,
        "n_null_draws": n_null,
        "methodology_note": (
            "STRING functional links; primary edge null draws random gene sets of the same cardinality "
            "as the high-burden union subset uniformly without replacement from network nodes (excluding "
            "fixed mechanism seeds) when the pool is large enough; otherwise falls back to degree-bin "
            "matched draws. Stress-distance null redraws random gene sets of the same size as "
            "PROGENy-filtered stress genes (not pathway-preserving)."
        ),
        "caveat": "Association network; not causal; induced subgraph is limited to posted identifiers.",
    }
    def _json_safe(o: object) -> object:
        if isinstance(o, float) and (math.isnan(o) or math.isinf(o)):
            return None
        if isinstance(o, dict):
            return {k: _json_safe(v) for k, v in o.items()}
        return o

    (out_dir / "exploratory_string_piezo1_bridge_summary.json").write_text(
        json.dumps(_json_safe(summ), indent=2), encoding="utf-8"
    )
    log(json.dumps(summ, indent=2))
