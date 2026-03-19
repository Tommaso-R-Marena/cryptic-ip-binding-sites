.PHONY: default install test lint format clean

default: text

install:
	pip install -r requirements.txt
	pip install pytest pytest-cov black flake8 mypy pandas matplotlib seaborn pre-commit

test:
	PYTHONPATH=. pytest tests/ --cov=pipeline

lint:
	flake8 . --count --select=E9,F63,F7,F82 --show-source --statistics
	mypy pipeline/

format:
	black .

clean:
	find . -type d -name "__pycache__" -exec rm -rf {} +
	find . -type d -name ".pytest_cache" -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
