"""Resource and memory controls for long-running analyses."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable


def set_memory_limit(max_gb: float) -> None:
    """Set soft memory cap (GB) on POSIX systems; no-op elsewhere."""
    try:
        import resource

        max_bytes = int(max_gb * 1024**3)
        resource.setrlimit(resource.RLIMIT_AS, (max_bytes, max_bytes))
    except Exception:
        return


def cleanup_files(paths: Iterable[Path]) -> None:
    """Delete intermediate files/directories if they exist."""
    for path in paths:
        if not path.exists():
            continue
        if path.is_dir():
            for child in path.glob("**/*"):
                if child.is_file():
                    child.unlink(missing_ok=True)
            for child in sorted(path.glob("**/*"), reverse=True):
                if child.is_dir():
                    child.rmdir()
            path.rmdir()
        else:
            path.unlink(missing_ok=True)
