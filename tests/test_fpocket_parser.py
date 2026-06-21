"""Tests for fpocket output parsing."""

from pathlib import Path

from cryptic_ip.analysis.fpocket_parser import FpocketParser


def test_fpocket_parser_extracts_residue_ids(tmp_path):
    info = tmp_path / "protein_info.txt"
    pockets_dir = tmp_path / "pockets"
    pockets_dir.mkdir()
    pocket_file = pockets_dir / "pocket1_atm.pdb"
    pocket_file.write_text(
        "ATOM    433   CB ARG A  58      -3.431  -1.275  -3.444  0.00  0.00           C 0\n"
        "ATOM    991   CA TYR A 128       6.185  13.322  -5.058  0.00  0.00           C 0\n",
        encoding="utf-8",
    )
    info.write_text(
        "Pocket 1 :\n"
        "Score : 0.42\n"
        "Volume : 120.0\n",
        encoding="utf-8",
    )

    frame = FpocketParser().parse_info_file(info)
    assert not frame.empty
    assert frame.iloc[0]["fpocket_residue_ids"] == "58,128"
