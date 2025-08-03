# Makefile for olmsted-cli

.PHONY: test clean clean-test clean-pyc clean-build help

help:
	@echo "Available commands:"
	@echo "  make test        - Run pytest"
	@echo "  make clean       - Remove all build, test, and Python artifacts"
	@echo "  make clean-test  - Remove test artifacts"
	@echo "  make clean-pyc   - Remove Python cache files"
	@echo "  make clean-build - Remove build artifacts"

test:
	pytest

clean: clean-build clean-pyc clean-test

clean-test:
	rm -rf _test_output/
	rm -rf .pytest_cache/
	rm -rf .coverage
	rm -rf htmlcov/
	rm -rf .tox/

clean-pyc:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -rf {} +

clean-build:
	rm -rf build/
	rm -rf dist/
	rm -rf .eggs/
	find . -name '*.egg-info' -exec rm -rf {} +
	find . -name '*.egg' -exec rm -f {} +