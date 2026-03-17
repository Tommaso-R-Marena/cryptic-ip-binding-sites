from pathlib import Path

import pandas as pd

from cryptic_ip.reproducibility import deterministic_sort_dataframe, load_yaml, validate_config


def test_default_config_validates_against_schema():
    config = load_yaml(Path("config/defaults/pipeline.yaml"))
    validate_config(config, Path("config/schemas/pipeline.schema.yaml"))


def test_deterministic_sort_dataframe_with_tie_breakers():
    frame = pd.DataFrame(
        [
            {"composite_score": 0.9, "uniprot_id": "B", "pocket_id": 2},
            {"composite_score": 0.9, "uniprot_id": "A", "pocket_id": 5},
            {"composite_score": 0.95, "uniprot_id": "Z", "pocket_id": 1},
        ]
    )

    sorted_frame = deterministic_sort_dataframe(
        frame,
        ["composite_score", "uniprot_id", "pocket_id"],
        ascending=[False, True, True],
    )

    assert list(sorted_frame["uniprot_id"]) == ["Z", "A", "B"]
