# Contributing to Cryptic IP Binding Sites

Thank you for your interest in contributing! This is an active research project.

## Ways to Contribute

### 1. Report Bugs

Found a bug? Please [open an issue](https://github.com/Tommaso-R-Marena/cryptic-ip-binding-sites/issues) with:
- Clear description of the problem
- Steps to reproduce
- Expected vs. actual behavior
- System information (OS, Python version)
- Relevant error messages or logs

### 2. Suggest Features

Have an idea? Open an issue with:
- Clear description of the feature
- Use cases and examples
- Why it would be valuable

### 3. Contribute Code

#### Development Setup

```bash
# Fork and clone
git clone https://github.com/YOUR-USERNAME/cryptic-ip-binding-sites.git
cd cryptic-ip-binding-sites

# Create development environment
conda env create -f environment.yml
conda activate cryptic-ip

# Install in editable mode with dev dependencies
pip install -e ".[dev]"

# Install pre-commit hooks
pre-commit install
```

#### Code Standards

- **Style**: Follow PEP 8, enforced by `black` and `flake8`
- **Type hints**: Use type annotations for all functions
- **Documentation**: Docstrings for all public functions (NumPy style)
- **Testing**: Write tests for new features
- **Commits**: Clear, descriptive commit messages

#### Testing

```bash
# Run full test suite
pytest tests/ -v

# Run with coverage
pytest tests/ --cov=cryptic_ip --cov-report=html

# Run specific test file
pytest tests/test_analysis.py -v

# Test notebooks
jupyter nbconvert --execute notebooks/*.ipynb
```

#### Pull Request Process

1. Create a feature branch: `git checkout -b feature/my-feature`
2. Make your changes
3. Run tests and linting: `pytest tests/ && black . && flake8`
4. Commit with clear messages
5. Push to your fork: `git push origin feature/my-feature`
6. Open a pull request

**PR Checklist:**
- [ ] Tests pass locally
- [ ] New tests added for new features
- [ ] Documentation updated
- [ ] Code follows style guidelines
- [ ] Commit messages are clear

### 4. Improve Documentation

Help improve docs by:
- Fixing typos or unclear sections
- Adding examples
- Improving API documentation
- Writing tutorials

### 5. Validate Results

Help validate the pipeline:
- Test on known IP-binding proteins
- Compare results with literature
- Suggest improvements to scoring

## Questions?

Feel free to open an issue or contact:
- **Email**: marena@cua.edu
- **GitHub**: [@Tommaso-R-Marena](https://github.com/Tommaso-R-Marena)

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
