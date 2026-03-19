import pytest
import os
import shutil
from pathlib import Path

# Add a marker to run tests with full environment if needed, otherwise skip missing tools
def has_tool(name):
    return shutil.which(name) is not None
    
# Basic sanity check that pytest is working
def test_pytest_framework():
    assert True
