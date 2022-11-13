.DEFAULT_GOAL := help

.PHONY:
	clean
	dist
	lint
	upload
	upload-test
	help

clean:  ## remove all build, coverage, Python, and test artifacts
	rm -f .coverage
	rm -f .installed.cfg
	rm -f coverage.xml
	rm -f pytest.xml
	rm -fr .eggs/
	rm -fr .mypy_cache/
	rm -fr .pytest_cache/
	rm -fr .tox/
	rm -fr bin/
	rm -fr build/
	rm -fr develop-eggs/
	rm -fr dist/
	rm -fr eggs/
	rm -fr htmlcov/
	rm -fr parts/
	find . -name '*.egg' -exec rm -fr {} +
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

dist: clean  ## generate the distributable files
	python3 -m pip install --upgrade build
	python3 -m build

lint:  ## fix code style and formatting
	isort ensembl_map/ tests/
	black ensembl_map/ tests/
	mypy ensembl_map/ tests/
	flake8 ensembl_map/ tests/

upload: dist  ## upload package to PyPI
	twine upload dist/*

upload-test: dist  ## upload package to TestPyPI
	twine upload --repository testpypi dist/*

help:  ## show this message and exit
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "\033[32m%-13s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)
