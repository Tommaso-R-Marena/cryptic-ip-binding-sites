"""Profiling helpers for timing and memory diagnostics."""

from __future__ import annotations

import functools
import threading
import time
from dataclasses import dataclass
from typing import Any, Callable, Dict, List


@dataclass
class TimingRecord:
    """Aggregated runtime stats for a profiled callable."""

    function_name: str
    call_count: int
    total_seconds: float
    avg_seconds: float
    max_seconds: float


class TimingRegistry:
    """Thread-safe storage for function-level timing observations."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._records: Dict[str, List[float]] = {}

    def add(self, function_name: str, elapsed_seconds: float) -> None:
        with self._lock:
            self._records.setdefault(function_name, []).append(elapsed_seconds)

    def summary(self) -> List[TimingRecord]:
        with self._lock:
            result: List[TimingRecord] = []
            for function_name, durations in self._records.items():
                total = sum(durations)
                count = len(durations)
                result.append(
                    TimingRecord(
                        function_name=function_name,
                        call_count=count,
                        total_seconds=total,
                        avg_seconds=total / count,
                        max_seconds=max(durations),
                    )
                )
            return sorted(result, key=lambda item: item.total_seconds, reverse=True)

    def clear(self) -> None:
        with self._lock:
            self._records.clear()


TIMING_REGISTRY = TimingRegistry()


def timed(stage_name: str | None = None) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
    """Decorator that records wall-clock runtime for function calls."""

    def _decorate(func: Callable[..., Any]) -> Callable[..., Any]:
        name = stage_name or func.__qualname__

        @functools.wraps(func)
        def _wrapped(*args: Any, **kwargs: Any) -> Any:
            started = time.perf_counter()
            try:
                return func(*args, **kwargs)
            finally:
                TIMING_REGISTRY.add(name, time.perf_counter() - started)

        return _wrapped

    return _decorate
