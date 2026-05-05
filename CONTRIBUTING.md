# Contributing to HydroRaVENS

## Reporting bugs and requesting features

Please open an issue on [GitHub](https://github.com/MNiMORPH/hydroRaVENS/issues).
Include a minimal reproducible example and your Python and package versions where relevant.

## Contributing code

1. Fork the repository and create a branch from `master`.
2. Install the package in editable mode with documentation dependencies:
   ```
   pip install -e ".[docs]"
   ```
3. Make your changes. Run the Cannon River example to confirm the model still produces a sensible result:
   ```
   cd examples/cannon_river
   python driver.py
   ```
4. If you changed the documentation, build it locally to check for errors:
   ```
   cd docs
   make html
   ```
5. Open a pull request against `master` with a clear description of what changed and why.

## Code style

- Follow [PEP 8](https://peps.python.org/pep-0008/). Use `pycodestyle --max-line-length=99` to check.
- Write [NumPy-style docstrings](https://numpydoc.readthedocs.io/en/latest/format.html) for all public methods.
- Intentional PEP 8 exceptions in this codebase: E741 (single-letter scientific variables such as `I`), E221 (alignment spaces in assignments).

## Contact

Andrew D. Wickert — [awickert@umn.edu](mailto:awickert@umn.edu) — [ORCID 0000-0002-9545-3365](https://orcid.org/0000-0002-9545-3365)
