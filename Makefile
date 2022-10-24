.DEFAULT_GOAL := help

.PHONY:
	clean
	install
	install-dev
	lint
	test
	upload
	venv
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

install: clean  ## install the package for production
	pip install -U pip
	pip install -U setuptools
	pip install -U .

install-dev: clean  ## install the package for development
	pip install -U pip
	pip install -U setuptools
	pip install -U -e .[dev]

lint:  ## fix code style and formatting
	isort --line-length 100 --profile black ensembl_map/ tests/ setup.py
	black -C --line-length 100 ensembl_map/ tests/ setup.py
	mypy --ignore-missing-imports ensembl_map/ tests/ setup.py
	flake8 --ignore E203,E501,W291,W503 ensembl_map/ tests/ setup.py

test:  ## run tests
	coverage run --source ensembl_map -m pytest --junitxml=pytest.xml
	coverage report -m
	coverage html
	coverage xml

test-cache:  # build cache required for tests
	python tests/make_cache.py

upload: clean  ## upload package to the GSC pypi server
	python setup.py sdist
	python setup.py install
	python setup.py bdist_wheel
	twine upload dist/*

help:  ## show this message and exit
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "\033[32m%-13s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)
