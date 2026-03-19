import pytest
import pandas as pd
from unittest.mock import patch, MagicMock
from screening.triage import triage_results

@patch("screening.triage.plt")
@patch("screening.triage.sns")
@patch("screening.triage.pd.read_csv")
@patch("screening.triage.Path.exists")
def test_triage_results(mock_exists, mock_read_csv, mock_sns, mock_plt):
    mock_exists.return_value = True
    
    # Mocking the master CSV dataframe
    mock_df = pd.DataFrame({
        'UniProtID': ['P12345', 'Q67890', 'A0A0B4J2F0'],
        'ManualReviewFlag': [True, False, True],
        'CompositeScore': [0.85, 0.40, 0.75],
        'MeanSASA': [15.0, 50.0, 20.0],
        'MeanPLDDT': [90.0, 40.0, 85.0]
    })
    
    mock_read_csv.return_value = mock_df
    
    # Patch DataFrame to_csv to avoid writing to disk during test
    with patch.object(pd.DataFrame, 'to_csv') as mock_to_csv:
        triage_results("test_organism")
        
    mock_read_csv.assert_called_once()
    
    # Should have saved the candidates CSV
    mock_to_csv.assert_called_once()
    
    # Should have plotted the score distribution and scatter plot
    assert mock_sns.histplot.called
    assert mock_sns.scatterplot.called
    assert mock_plt.savefig.call_count == 2
