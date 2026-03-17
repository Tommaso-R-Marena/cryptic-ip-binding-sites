"""Input validation utilities with optional Pydantic-backed models."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional

from ..errors import UnsupportedFormatError, ValidationError

try:
    from pydantic import (
        BaseModel,
        ConfigDict,
        Field,
        ValidationError as PydanticValidationError,
        field_validator,
    )

    PYDANTIC_AVAILABLE = True
except Exception:  # pragma: no cover - environment fallback
    BaseModel = object  # type: ignore[assignment]
    ConfigDict = dict  # type: ignore[assignment]
    Field = lambda default=None, **kwargs: default  # type: ignore[assignment]
    PydanticValidationError = ValueError
    PYDANTIC_AVAILABLE = False

    def field_validator(*_args, **_kwargs):
        def _decorator(func):
            return func

        return _decorator


SUPPORTED_STRUCTURE_EXTENSIONS = {".pdb", ".cif", ".ent"}
SUPPORTED_EXPORT_FORMATS = {"csv", "json", "hdf5"}


class BatchDownloadConfig(BaseModel):
    if PYDANTIC_AVAILABLE:
        model_config = ConfigDict(extra="forbid")
    requests_per_second: float = Field(default=4.0)
    max_retries: int = Field(default=5)
    backoff_seconds: float = Field(default=1.5)
    timeout_seconds: float = Field(default=30.0)


class ParallelProcessorConfig(BaseModel):
    if PYDANTIC_AVAILABLE:
        model_config = ConfigDict(extra="forbid")
    workers: int = Field(default=1)
    chunk_size: int = Field(default=100)
    checkpoint_every: int = Field(default=1)


class ExportConfig(BaseModel):
    if PYDANTIC_AVAILABLE:
        model_config = ConfigDict(extra="forbid")
    export_format: str

    @field_validator("export_format")
    @classmethod
    def _valid_format(cls, value: str) -> str:
        normalized = value.lower()
        if normalized not in SUPPORTED_EXPORT_FORMATS:
            raise ValueError(
                f"Unsupported export format '{value}'. Use one of: {', '.join(sorted(SUPPORTED_EXPORT_FORMATS))}."
            )
        return normalized


def validate_structure_file(path: Path | str, allowed: Optional[Iterable[str]] = None) -> Path:
    file_path = Path(path)
    if not file_path.exists():
        raise ValidationError(f"Input file not found: {file_path}")

    allowed_suffixes = {suffix.lower() for suffix in (allowed or SUPPORTED_STRUCTURE_EXTENSIONS)}
    suffix = file_path.suffix.lower()
    if suffix not in allowed_suffixes:
        raise UnsupportedFormatError(file_path, ", ".join(sorted(allowed_suffixes)))

    return file_path


def _fallback_validate(model_class: type[BaseModel], payload: dict):
    obj = model_class()
    for key, value in payload.items():
        setattr(obj, key, value)

    if isinstance(obj, BatchDownloadConfig):
        if not (0.1 <= obj.requests_per_second <= 50.0):
            raise ValidationError("requests_per_second must be between 0.1 and 50.0")
        if not (1 <= obj.max_retries <= 20):
            raise ValidationError("max_retries must be between 1 and 20")
        if not (0.1 <= obj.backoff_seconds <= 30.0):
            raise ValidationError("backoff_seconds must be between 0.1 and 30.0")
        if not (1.0 <= obj.timeout_seconds <= 600.0):
            raise ValidationError("timeout_seconds must be between 1.0 and 600.0")

    if isinstance(obj, ParallelProcessorConfig):
        if not (1 <= obj.workers <= 512):
            raise ValidationError("workers must be between 1 and 512")
        if not (1 <= obj.chunk_size <= 10000):
            raise ValidationError("chunk_size must be between 1 and 10000")
        if not (1 <= obj.checkpoint_every <= 5000):
            raise ValidationError("checkpoint_every must be between 1 and 5000")

    if isinstance(obj, ExportConfig):
        obj.export_format = obj._valid_format(obj.export_format)

    return obj


def parse_or_raise(model_class: type[BaseModel], payload: dict, message: str):
    try:
        if PYDANTIC_AVAILABLE:
            return model_class(**payload)
        return _fallback_validate(model_class, payload)
    except (PydanticValidationError, ValueError, TypeError, ValidationError) as exc:
        raise ValidationError(f"{message}: {exc}") from exc
