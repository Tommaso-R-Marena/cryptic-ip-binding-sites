"""
Database module for proteome data management.
"""

from .downloader import ProteomeDownloader
from .manager import ProteomeManager

__all__ = ["ProteomeDownloader", "ProteomeManager"]
