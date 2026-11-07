# Contributing to twas-xcal

Thank you for your interest. Bug reports, questions and pull requests are welcome.

## Reporting a problem
Open an [issue](https://github.com/mqar96/twas-xcal/issues) with:
- what you ran (command or code snippet),
- the full error message or unexpected output,
- your Python version and operating system.

## Development setup
```bash
git clone https://github.com/mqar96/twas-xcal.git
cd twas-xcal
pip install -e ".[test]"
python -m pytest tests
```

## Pull requests
- Keep changes focused; one topic per pull request.
- Add or update tests for any change in statistical behaviour.
- Make sure `python -m pytest tests` passes.
- For methodological changes, describe the expected effect on type-I error and power.

## Contact
Muhammad Qudduse Ahmad Rao — qudduse1234@gmail.com — ORCID [0009-0001-8731-7080](https://orcid.org/0009-0001-8731-7080)
