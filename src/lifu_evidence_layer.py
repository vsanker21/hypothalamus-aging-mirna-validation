"""Write LIFU / tFUS evidence layer (kept separate from MOESM statistics)."""
from pathlib import Path

from config import OUTPUT_DIR

TEXT = """Evidence layer: LIFU / neuromodulation (separate from MOESM statistics)
-----------------------------------------------------------------------------
This file documents human and translational safety/physics rationale for
transcranial low-intensity focused ultrasound (LIFU / tFUS). It is NOT merged
statistically with Zhang et al. (2017) mouse supplementary tables in this
repository. Treat as parallel evidence streams for grant / IRB narratives.

1) Human wearable LIFU neuromodulation (precedent hardware class)
   - Bawiec et al., Journal of Ultrasound in Medicine 2025 (doi:10.1002/jum.16600):
     wearable steerable transcranial LIFU; volunteer study (n=20) with SWI;
     hydrophone/simulation agreement; skull/sinus steering discussion.

2) Systematic safety context for human tFUS neuromodulation
   - Sarica et al., Brain Stimulation 2022: systematic review of human tFUS.
   - ITRUSST consensus (Aubry / Ter Haar et al., arXiv:2311.05359, 2023):
     biophysical safety limits for transcranial ultrasonic stimulation.

3) Regulatory / device-class precedents (non-exhaustive)
   - MR-guided HIFU approvals (e.g. essential tremor) establish transcranial
     focusing precedent at ablative intensities; neuromodulation LIFU operates
     far below tissue-lesioning thresholds (different risk class).
   - Map milestones early to research-device vs clinical regulatory categories
     with institutional regulatory counsel before human-subjects language appears
     in proposals (see data/provenance/TRANSLATION_REGULATORY_FRAMING.txt).

4) What would be required to merge LIFU with MOESM statistically?
   - New experiments: defined LIFU dose to MBH/htNSC-relevant coordinates with sham
     controls, paired molecular readouts (e.g. exosome miRNA, NSC markers).

5) Acoustic reporting checklist (future work)
   - MI, ISPTA, frequency, PRF, duty cycle, focal depth, skull correction,
     thermometry if applicable (ITRUSST).

6) Narrative bridge (hypothesis planning only; not statistical evidence)
   - For Science Advances–style framing of hypothalamic stem/niche biology + aging omics +
     LIFU as a candidate neuromodulatory interface, see:
       data/provenance/FRAMING_NSC_HYPOTHALAMUS_LIFU_AGING_SCIENCE_ADVANCES.txt
   - That document lists figure panel titles, a claims ladder, and pre-experimental
     predictions/falsifiers. It does not authorize merging this LIFU lane with omics tests.

Cite this file alongside repo outputs with explicit data-layer attribution.
"""


def write_lifu_layer() -> Path:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    p = OUTPUT_DIR / "LIFU_evidence_layer.txt"
    p.write_text(TEXT, encoding="utf-8")
    return p
