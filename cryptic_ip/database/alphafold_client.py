"""AlphaFold Database API client for fetching protein structures."""

import logging
import re
import requests
from pathlib import Path
from typing import Dict, Optional, List

from .async_fetch import FetchJob, fetch_all, fetch_alphafold_models

logger = logging.getLogger(__name__)


class AlphaFoldClient:
    """Client for AlphaFold Protein Structure Database.

    Model files are located through the prediction API
    (https://alphafold.ebi.ac.uk/api-docs) and fetched by
    :mod:`cryptic_ip.database.async_fetch`. Whole proteomes should come from the
    per-proteome archives on the FTP site
    (https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/), as
    ``scripts/proteome_screen.py`` and the proteome-screen workflow do; the FTP
    ``latest`` directory holds only those archives, not individual models.
    """
    
    API_BASE = "https://alphafold.ebi.ac.uk/api"
    FTP_BASE = "https://ftp.ebi.ac.uk/pub/databases/alphafold/latest"
    
    def __init__(self, cache_dir: Optional[Path] = None):
        """Initialize AlphaFold client.
        
        Args:
            cache_dir: Directory to cache downloaded structures
        """
        self.cache_dir = cache_dir or Path.home() / ".cryptic_ip" / "alphafold_cache"
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'CrypticIP/1.0 (https://github.com/Tommaso-R-Marena/cryptic-ip-binding-sites)'
        })
    
    def _cached_versions(self, uniprot_id: str) -> List[Path]:
        """Cached models for ``uniprot_id``, newest version first.

        Sorted by the numeric version, not the file name: lexically "v10"
        sorts before "v4".
        """
        def version_of(path: Path) -> int:
            match = re.search(r"_v(\d+)\.pdb$", path.name)
            return int(match.group(1)) if match else -1

        found = self.cache_dir.glob(f"AF-{uniprot_id}-F1-model_v*.pdb")
        return sorted(found, key=version_of, reverse=True)

    def fetch_structure(self, uniprot_id: str, version: Optional[int] = None) -> Path:
        """Download the AlphaFold model for a UniProt accession.

        The current model's URL comes from the prediction API, so a newer
        release is fetched rather than a stale cached one being returned in its
        place. The download is retried on transient failures, checked to be a
        complete PDB file, and written atomically.

        Only when the API cannot be reached is the newest cached model returned,
        with a warning.

        Args:
            uniprot_id: UniProt accession.
            version: Fetch this model version specifically instead of the
                current one.

        Raises:
            ValueError: The accession has no AlphaFold model.
            ConnectionError: The model could not be retrieved and none is cached.
        """
        uniprot_id = uniprot_id.strip().upper()
        if version is not None:
            filename = f"AF-{uniprot_id}-F1-model_v{int(version)}.pdb"
            (result,) = fetch_all(
                [FetchJob(uniprot_id, f"https://alphafold.ebi.ac.uk/files/{filename}", self.cache_dir / filename)],
                concurrency=1,
            )
        else:
            result = fetch_alphafold_models([uniprot_id], self.cache_dir, concurrency=1)[uniprot_id]

        if result.ok:
            if not result.from_cache:
                logger.info("Downloaded %s to %s", uniprot_id, result.path)
            return Path(result.path)
        if result.not_found:
            raise ValueError(
                f"UniProt ID {uniprot_id} not found in AlphaFold Database "
                f"(https://alphafold.ebi.ac.uk/entry/{uniprot_id}): {result.error}"
            )
        cached = self._cached_versions(uniprot_id)
        if cached and version is None:
            logger.warning(
                "AlphaFold unreachable (%s); using cached %s, which may not be the current release",
                result.error,
                cached[0].name,
            )
            return cached[0]
        raise ConnectionError(f"Could not retrieve AlphaFold model for {uniprot_id}: {result.error}")

    def get_metadata(self, uniprot_id: str) -> Dict:
        """Fetch AlphaFold prediction metadata.
        
        Args:
            uniprot_id: UniProt accession
            
        Returns:
            Dictionary with prediction metadata
        """
        url = f"{self.API_BASE}/prediction/{uniprot_id}"
        
        try:
            response = self.session.get(url, timeout=10)
            response.raise_for_status()
            data = response.json()
            
            if not data:
                raise ValueError(f"No AlphaFold prediction for {uniprot_id}")
            
            # Extract first entry (usually only one)
            entry = data[0] if isinstance(data, list) else data
            
            return {
                'uniprot_id': entry['uniprotAccession'],
                'entry_id': entry.get('entryId', f"AF-{uniprot_id}-F1"),
                'gene': entry.get('gene', ''),
                'organism': entry.get('organismScientificName', ''),
                'sequence_length': entry.get('uniprotSequenceLength', 0),
                'global_metric_value': entry.get('globalMetricValue', 0.0),
                'model_version': entry.get('latestVersion', 4),
                'model_date': entry.get('modelCreatedDate', ''),
            }
            
        except requests.HTTPError as e:
            if e.response.status_code == 404:
                raise ValueError(f"No metadata for {uniprot_id}")
            raise
    
    def fetch_batch(
        self, uniprot_ids: List[str], delay: float = 0.0, concurrency: int = 16
    ) -> Dict[str, Optional[Path]]:
        """Download many models concurrently.

        Args:
            uniprot_ids: UniProt accessions.
            delay: Minimum interval between requests to one host, in seconds.
            concurrency: Downloads in flight at once.

        Returns:
            Mapping accession -> model path, or ``None`` when unavailable.
        """
        results = fetch_alphafold_models(
            uniprot_ids, self.cache_dir, concurrency=concurrency, min_interval=delay,
            manifest=self.cache_dir / "download_manifest.json",
        )
        return {key: Path(r.path) if r.ok else None for key, r in results.items()}

    def fetch_proteome(self, proteome_id: str, output_dir: Path) -> int:
        """Download entire AlphaFold proteome.
        
        Note: Large proteomes (>10 GB). Download from FTP is recommended:
        https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/
        
        Args:
            proteome_id: UniProt proteome ID (e.g., 'UP000002311' for yeast)
            output_dir: Directory to save structures
            
        Returns:
            Number of structures downloaded
        """
        logger.warning(
            f"Downloading full proteome {proteome_id}. "
            f"This may take hours. Consider using FTP bulk download:"
            f"\n  wget https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/{proteome_id}_*.tar"
        )
        
        # For now, direct users to FTP
        raise NotImplementedError(
            "Bulk proteome download not implemented via API. "
            "Please use FTP download as documented in README.md"
        )


if __name__ == "__main__":
    # Example usage
    logging.basicConfig(level=logging.INFO)
    
    client = AlphaFoldClient()
    
    print("Testing AlphaFold client...\n")
    
    # Test 1: Get metadata
    print("1. Fetching ADAR2 metadata...")
    try:
        metadata = client.get_metadata('P78563')
        print(f"   Gene: {metadata['gene']}")
        print(f"   Organism: {metadata['organism']}")
        print(f"   Length: {metadata['sequence_length']} aa")
        print(f"   ✓ Metadata OK\n")
    except Exception as e:
        print(f"   ✗ Metadata failed: {e}\n")
    
    # Test 2: Download structure
    print("2. Downloading ADAR2 structure...")
    try:
        adar2_path = client.fetch_structure('P78563')
        print(f"   File: {adar2_path}")
        print(f"   Size: {adar2_path.stat().st_size / 1024:.1f} KB")
        print(f"   ✓ Download OK\n")
    except Exception as e:
        print(f"   ✗ Download failed: {e}\n")
    
    print("Testing complete!")
