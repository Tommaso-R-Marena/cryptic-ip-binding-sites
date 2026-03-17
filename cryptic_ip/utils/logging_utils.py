"""Structured logging helpers for pipeline and batch execution."""

from __future__ import annotations

import json
import logging
import logging.handlers
import os
from pathlib import Path
from typing import Any, Dict, Optional


class JsonFormatter(logging.Formatter):
    """Minimal JSON formatter suitable for file aggregation and journald/syslog ingestion."""

    def format(self, record: logging.LogRecord) -> str:
        payload: Dict[str, Any] = {
            "timestamp": self.formatTime(record, self.datefmt),
            "level": record.levelname,
            "logger": record.name,
            "message": record.getMessage(),
        }
        if hasattr(record, "context"):
            payload["context"] = getattr(record, "context")
        if record.exc_info:
            payload["exception"] = self.formatException(record.exc_info)
        return json.dumps(payload, default=str)


def configure_logging(
    log_dir: Path | str,
    module_name: str,
    level: int = logging.INFO,
    enable_syslog: bool = False,
) -> logging.Logger:
    log_path = Path(log_dir)
    log_path.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger(module_name)
    logger.setLevel(level)
    logger.propagate = False

    if logger.handlers:
        return logger

    formatter = JsonFormatter()

    file_handler = logging.FileHandler(log_path / f"{module_name}.log")
    file_handler.setLevel(level)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    aggregate_handler = logging.FileHandler(log_path / "batch_jobs.log")
    aggregate_handler.setLevel(level)
    aggregate_handler.setFormatter(formatter)
    logger.addHandler(aggregate_handler)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(level)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    if enable_syslog:
        _attach_system_handler(logger, formatter, level)

    return logger


def _attach_system_handler(logger: logging.Logger, formatter: logging.Formatter, level: int) -> None:
    """Attach journald/syslog handlers when available in runtime environment."""
    journald_socket = "/run/systemd/journal/socket"
    if os.path.exists(journald_socket):
        handler = logging.handlers.SysLogHandler(address=journald_socket)
    else:
        handler = logging.handlers.SysLogHandler(address=("localhost", 514))

    handler.setLevel(level)
    handler.setFormatter(formatter)
    logger.addHandler(handler)


def log_with_context(logger: logging.Logger, level: int, message: str, **context: Any) -> None:
    logger.log(level, message, extra={"context": context})
