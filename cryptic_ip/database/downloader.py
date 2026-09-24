"""
Download and organize AlphaFold proteome structures.
"""

import re
import tarfile
from pathlib import Path
from typing import Dict

from tqdm import tqdm

from .async_fetch import FetchJob, download_large_file, fetch_all


def newest_archive(listing: str, proteome_id: str) -> str:
    """The newest ``{proteome}_{taxon}_{MNEMONIC}_v{n}.tar`` named in a listing."""
    found = set(re.findall(rf"({re.escape(proteome_id)}_\d+_[A-Z0-9]+_v\d+\.tar)", listing))
    if not found:
        raise RuntimeError(f"No archive for {proteome_id} in the listing")
    return max(found, key=lambda name: int(re.search(r"_v(\d+)\.tar$", name).group(1)))


def extract_models(tar_path: Path, dest: Path, *, pattern: str = ".pdb.gz") -> int:
    """Extract the archive members ending in ``pattern`` into ``dest``.

    Members are extracted by base name only, so an archive entry with a path
    component cannot write outside ``dest``; the mmCIF copies AlphaFold ships
    alongside are skipped unless asked for.
    """
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    count = 0
    with tarfile.open(tar_path, "r:*") as archive:
        for member in archive:
            name = Path(member.name).name
            if not member.isfile() or not name.endswith(pattern):
                continue
            source = archive.extractfile(member)
            if source is None:
                continue
            target = dest / name
            tmp = target.with_name(f".{name}.part")
            with open(tmp, "wb") as handle:
                while True:
                    block = source.read(1 << 20)
                    if not block:
                        break
                    handle.write(block)
            tmp.replace(target)
            count += 1
    return count


class ProteomeDownloader:
    """
    Download AlphaFold proteome structures.
    
    Supported organisms:
    - Saccharomyces cerevisiae (yeast): UP000002311
    - Homo sapiens (human): UP000005640
    - Dictyostelium discoideum: UP000002195
    """
    
    PROTEOMES = {
        'yeast': {
            'uniprot_id': 'UP000002311',
            'organism': 'Saccharomyces cerevisiae',
            'taxon': '559292',
            'proteins': 6049,
            'size_gb': 15
        },
        'human': {
            'uniprot_id': 'UP000005640',
            'organism': 'Homo sapiens',
            'taxon': '9606',
            'proteins': 23391,
            'size_gb': 50
        },
        'dictyostelium': {
            'uniprot_id': 'UP000002195',
            'organism': 'Dictyostelium discoideum',
            'taxon': '44689',
            'proteins': 12622,
            'size_gb': 30
        }
    }
    
    ALPHAFOLD_BASE = "https://ftp.ebi.ac.uk/pub/databases/alphafold/latest"
    
    def __init__(self, data_dir: str = "data/structures"):
        """
        Initialize downloader.
        
        Args:
            data_dir: Base directory for proteome structures
        """
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)
    
    def download_proteome(self, organism: str, force: bool = False) -> Path:
        """
        Download complete proteome from AlphaFold.
        
        Args:
            organism: Organism key ('yeast', 'human', 'dictyostelium')
            force: Re-download even if exists
            
        Returns:
            Path to proteome directory
        """
        if organism not in self.PROTEOMES:
            raise ValueError(f"Unknown organism: {organism}")

        info = self.PROTEOMES[organism]
        proteome_dir = self.data_dir / organism
        if proteome_dir.exists() and any(proteome_dir.glob("AF-*.pdb*")) and not force:
            print(f"{organism} proteome already exists at {proteome_dir}")
            return proteome_dir
        proteome_dir.mkdir(parents=True, exist_ok=True)

        # The archive name carries the release version and a mnemonic that is
        # not the organism key (Dictyostelium is "DICDI"), so it is read from
        # the directory listing. Building it from the key - as this method once
        # did, with a fixed "_v4" - fails for every proteome once v4 is gone.
        url = self.resolve_archive_url(organism)
        tar_path = self.data_dir / url.rsplit("/", 1)[-1]
        print(f"\nDownloading {organism} proteome: {url}")

        with tqdm(unit="B", unit_scale=True, desc=tar_path.name) as pbar:
            def report(done: int, total: int) -> None:
                pbar.total = total or None
                pbar.n = done
                pbar.refresh()

            result = download_large_file(url, tar_path, progress=report)
        if not result.ok:
            raise RuntimeError(f"Download of {url} failed: {result.error}")

        print(f"\nExtracting PDB-format models to {proteome_dir}...")
        try:
            n = extract_models(tar_path, proteome_dir, pattern=".pdb.gz")
        finally:
            tar_path.unlink(missing_ok=True)
        print(f"Extracted {n} models: {proteome_dir}")
        return proteome_dir

    def resolve_archive_url(self, organism: str) -> str:
        """URL of the newest archive for ``organism`` in the ``latest`` listing."""
        proteome_id = self.PROTEOMES[organism]["uniprot_id"]
        (listing,) = fetch_all([FetchJob("listing", f"{self.ALPHAFOLD_BASE}/")], concurrency=1)
        if not listing.ok:
            raise RuntimeError(f"Could not list {self.ALPHAFOLD_BASE}: {listing.error}")
        return f"{self.ALPHAFOLD_BASE}/{newest_archive(listing.payload.decode(errors='replace'), proteome_id)}"

    def get_info(self, organism: str) -> Dict:
        """
        Get information about a proteome.
        
        Args:
            organism: Organism key
            
        Returns:
            Proteome information dictionary
        """
        if organism not in self.PROTEOMES:
            raise ValueError(f"Unknown organism: {organism}")
        return self.PROTEOMES[organism]
    
    def list_available(self) -> Dict:
        """
        List all available proteomes.
        
        Returns:
            Dictionary of proteome information
        """
        return self.PROTEOMES
