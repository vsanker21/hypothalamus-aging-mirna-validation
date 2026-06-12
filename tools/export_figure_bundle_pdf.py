"""
Export sa_bundle PNG panels to 300 dpi PDF (one PDF per panel + combined supplement).

Run from feasibility_study/:
  python tools/export_figure_bundle_pdf.py --outputs-dir outputs
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BUNDLE = PROJECT_ROOT / "outputs" / "figures" / "sa_bundle"
TARGET_DPI = 300


def _png_to_pdf_page(png_path: Path, pdf_path: Path, dpi: int = TARGET_DPI) -> bool:
    if not png_path.is_file():
        return False
    img = Image.open(png_path)
    if img.mode not in ("RGB", "L"):
        img = img.convert("RGB")
    w_px, h_px = img.size
    fig_w = w_px / dpi
    fig_h = h_px / dpi
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)
    ax.imshow(img, aspect="auto")
    ax.axis("off")
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf_path, format="pdf", dpi=dpi, bbox_inches="tight", pad_inches=0)
    plt.close(fig)
    return True


def export_bundle(bundle_dir: Path, out_pdf_dir: Path, dpi: int = TARGET_DPI) -> tuple[int, Path | None]:
    bundle_dir = bundle_dir.resolve()
    out_pdf_dir = out_pdf_dir.resolve()
    out_pdf_dir.mkdir(parents=True, exist_ok=True)
    pngs = sorted(bundle_dir.glob("*.png"))
    if not pngs:
        return 0, None

    ok = 0
    for png in pngs:
        pdf_path = out_pdf_dir / f"{png.stem}.pdf"
        if _png_to_pdf_page(png, pdf_path, dpi=dpi):
            ok += 1

    combined = out_pdf_dir / "SUPPLEMENTARY_FIGURES_COMBINED.pdf"
    with PdfPages(combined) as pdf:
        for png in pngs:
            img = Image.open(png)
            if img.mode not in ("RGB", "L"):
                img = img.convert("RGB")
            w_px, h_px = img.size
            fig, ax = plt.subplots(figsize=(w_px / dpi, h_px / dpi), dpi=dpi)
            ax.imshow(img, aspect="auto")
            ax.axis("off")
            fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
            pdf.savefig(fig, bbox_inches="tight", pad_inches=0)
            plt.close(fig)

    manifest = out_pdf_dir / "FIGURE_PDF_MANIFEST.txt"
    manifest.write_text(
        "\n".join(
            [
                f"Figure PDF export — {len(pngs)} source PNG(s) at {dpi} dpi embedding",
                f"Bundle source: {bundle_dir}",
                "",
                *[f"  {p.name}" for p in sorted(out_pdf_dir.glob("*.pdf"))],
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return ok, combined


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outputs-dir", type=Path, default=PROJECT_ROOT / "outputs")
    ap.add_argument("--dpi", type=int, default=TARGET_DPI)
    args = ap.parse_args()
    bundle = args.outputs_dir / "figures" / "sa_bundle"
    out_dir = args.outputs_dir / "figures" / "sa_bundle_pdf"
    n, combined = export_bundle(bundle, out_dir, dpi=args.dpi)
    if n == 0:
        print(f"No PNG panels in {bundle}; run figures/build_sa_figure_bundle.py first.")
        return 1
    print(f"Exported {n} panel PDF(s) -> {out_dir}")
    if combined:
        print(f"Combined supplement: {combined}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
