# === FILE: docs/changelog.md ===
# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2023-10-27

### Added
- Complete Sphinx documentation system with `furo` theme.
- Google-style docstrings for all pipeline modules.
- Type hints for all public functions.
- Automated CI/CD workflow for documentation via GitHub Actions.
- Narrative documentation including scientific justification and workflow diagrams.
- API reference autogeneration using `sphinx-autoapi`.

### Changed
- Standardized logging and error handling across modules.
- Updated `requirements.txt` with documentation dependencies.

### Fixed
- Improved resilience of API calls in `06_analysis.py`.
