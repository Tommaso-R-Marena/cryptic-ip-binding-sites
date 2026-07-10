"""Colab runtime helpers: conda PATH and Python executable discovery."""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Optional

CONDA_PREFIXES: tuple[Path, ...] = (
    Path("/usr/local/miniforge3"),
    Path("/opt/conda"),
    Path("/usr/share/miniconda"),
    Path.home() / "miniforge3",
    Path.home() / "miniconda3",
)


def is_python_executable(path: str | Path) -> bool:
    """Return True when ``path`` points to a Python interpreter binary."""
    return Path(path).name.startswith("python")


def _resolve_python_exe(bin_dir: Path) -> Optional[Path]:
    """Pick the best Python executable inside a conda prefix."""
    for name in ("python", "python3"):
        candidate = bin_dir / name
        if candidate.exists():
            return candidate
    versioned = sorted(p for p in bin_dir.glob("python3.*") if p.is_file())
    return versioned[0] if versioned else None


def find_conda_prefix() -> Optional[Path]:
    """Return the first conda prefix that has fpocket installed."""
    for prefix in CONDA_PREFIXES:
        if (prefix / "bin" / "fpocket").exists():
            return prefix
    return None


def bootstrap_colab_runtime() -> str:
    """Put conda structural-biology tools on PATH; return pipeline Python executable."""
    prefix = find_conda_prefix()
    if prefix is None:
        return sys.executable

    bin_dir = prefix / "bin"
    path_entries = [str(bin_dir)]
    existing = os.environ.get("PATH", "")
    if str(bin_dir) not in existing.split(os.pathsep):
        os.environ["PATH"] = os.pathsep.join(path_entries + [existing])

    python_exe = _resolve_python_exe(bin_dir)
    if python_exe is not None:
        return str(python_exe)
    return sys.executable


def reexec_if_needed() -> None:
    """Re-launch the current script under conda Python when Colab uses system Python."""
    target = bootstrap_colab_runtime()
    if os.path.normpath(target) != os.path.normpath(sys.executable):
        os.execv(target, [target, *sys.argv])


def write_runtime_marker(repo_root: Path) -> Path:
    """Persist resolved python path for notebooks (written by colab_install.sh)."""
    marker_dir = repo_root / ".colab"
    marker_dir.mkdir(parents=True, exist_ok=True)
    python_path = marker_dir / "python_path"
    python_path.write_text(bootstrap_colab_runtime() + "\n", encoding="utf-8")
    return python_path


def read_runtime_python(repo_root: Path) -> str:
    """Read pipeline Python from marker file or discover conda runtime."""
    marker = repo_root / ".colab" / "python_path"
    if marker.exists():
        text = marker.read_text(encoding="utf-8").strip()
        if text and Path(text).exists():
            return text
    return bootstrap_colab_runtime()


if __name__ == "__main__":
    ROOT = Path(__file__).resolve().parents[1]
    marker = write_runtime_marker(ROOT)
    print(marker.read_text().strip())
