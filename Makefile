# Makefile for olmsted-cli

.PHONY: test schemas clean clean-test clean-pyc clean-build help

SCHEMA_DIR = olmsted_cli/schemas
YAML_SCHEMA = $(SCHEMA_DIR)/olmsted-schema.yaml
JSON_SCHEMA = $(SCHEMA_DIR)/olmsted-schema.json

help:
	@echo "Available commands:"
	@echo "  make test        - Run pytest"
	@echo "  make schemas     - Regenerate olmsted-schema.json from olmsted-schema.yaml"
	@echo "  make clean       - Remove all build, test, and Python artifacts"
	@echo "  make clean-test  - Remove test artifacts"
	@echo "  make clean-pyc   - Remove Python cache files"
	@echo "  make clean-build - Remove build artifacts"

test:
	pytest

# Regenerate the published JSON schema from the YAML source.
# Run this after editing olmsted-schema.yaml.
schemas: $(JSON_SCHEMA)

$(JSON_SCHEMA): $(YAML_SCHEMA)
	python3 -c "import yaml, json; data = yaml.safe_load(open('$<')); open('$@', 'w').write(json.dumps(data, indent=2) + '\n'); print('Generated $@')"

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