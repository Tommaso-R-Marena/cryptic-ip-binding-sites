"""Custom exception hierarchy for cryptic IP pipeline failures."""

from __future__ import annotations

from pathlib import Path
from typing import Optional


class CrypticIPError(Exception):
    """Base class for all pipeline-specific errors."""


class ValidationError(CrypticIPError):
    """Raised when user or file input is invalid."""


class NetworkRetryError(CrypticIPError):
    """Raised when a network request exhausts retry budget."""


class OperationTimeoutError(CrypticIPError):
    """Raised when an operation exceeds a configured timeout."""


class CacheOperationError(CrypticIPError):
    """Raised when cache/database operations fail permanently."""


class BatchItemProcessingError(CrypticIPError):
    """Raised when one item in a batch fails analysis."""

    def __init__(self, uniprot_id: str, reason: str):
        self.uniprot_id = uniprot_id
        self.reason = reason
        super().__init__(f"Failed to process {uniprot_id}: {reason}")


class UnsupportedFormatError(ValidationError):
    """Raised when a provided file does not match supported formats."""

    def __init__(self, path: Path, allowed: str):
        super().__init__(f"Unsupported file format for '{path.name}'. Allowed formats: {allowed}.")


class RecoveryStateError(CrypticIPError):
    """Raised when checkpoint/recovery state cannot be parsed safely."""


class ResourceLimitError(CrypticIPError):
    """Raised when runtime resources are insufficient for configured execution."""


class UserFacingError(CrypticIPError):
    """Exception with contextual, plain-language remediation."""

    def __init__(
        self,
        message: str,
        *,
        suggestion: Optional[str] = None,
        context: Optional[str] = None,
        docs_url: Optional[str] = None,
    ) -> None:
        self.suggestion = suggestion
        self.context = context
        self.docs_url = docs_url

        parts = [message]
        if context:
            parts.append(f"Context: {context}")
        if suggestion:
            parts.append(f"Suggested fix: {suggestion}")
        if docs_url:
            parts.append(f"Documentation: {docs_url}")
        super().__init__(" | ".join(parts))
