"""Write outputs/VERSION_SNAPSHOT.txt with Python + key R package versions (reproducibility)."""
from __future__ import annotations

import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    out = root / "outputs" / "VERSION_SNAPSHOT.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"UTC timestamp: {datetime.now(timezone.utc).isoformat()}",
        "",
        "=== Python ===",
        sys.version.replace("\n", " "),
        "",
        "=== pip freeze (project env) ===",
    ]
    try:
        proc = subprocess.run(
            [sys.executable, "-m", "pip", "freeze"],
            capture_output=True,
            text=True,
            timeout=120,
            encoding="utf-8",
            errors="replace",
        )
        lines.append(proc.stdout.strip() or proc.stderr.strip() or "(empty)")
    except Exception as e:
        lines.append(f"(pip freeze failed: {e})")

    lines.extend(["", "=== R ==="])
    try:
        rv = subprocess.run(
            ["Rscript", "-e", "cat(R.version.string, '\\n')"],
            capture_output=True,
            text=True,
            timeout=60,
            encoding="utf-8",
            errors="replace",
            cwd=str(root),
        )
        lines.append((rv.stdout or rv.stderr or "").strip())
    except Exception as e:
        lines.append(f"(Rscript not available: {e})")

    lines.extend(["", "=== R packages (Seurat edgeR limma fgsea) ==="])
    rcode = """
    for (p in c('Seurat','edgeR','limma','fgsea','BiocManager')) {
      if (requireNamespace(p, quietly=TRUE))
        message(p, ': ', as.character(packageVersion(p)))
      else
        message(p, ': NOT_INSTALLED')
    }
    """
    try:
        rv2 = subprocess.run(
            ["Rscript", "-e", rcode],
            capture_output=True,
            text=True,
            timeout=120,
            encoding="utf-8",
            errors="replace",
            cwd=str(root),
        )
        lines.append((rv2.stdout + rv2.stderr).strip())
    except Exception as e:
        lines.append(f"(R package probe failed: {e})")

    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(out.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
