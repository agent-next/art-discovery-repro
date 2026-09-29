.PHONY: setup check lint shellcheck test fmt

# Use the project venv when it exists (make setup / uv sync); fall back to PATH (CI).
BIN := $(if $(wildcard .venv/bin/ruff),.venv/bin/,)
COV_FLOOR := 94

setup:
	python3 -m venv .venv
	.venv/bin/pip install -e '.[dev]'

check: lint shellcheck test

lint:
	$(BIN)ruff check src tests pipeline benchmark

shellcheck:
	@if command -v shellcheck >/dev/null; then \
		shellcheck -S warning pipeline/*/*.sh; \
	else \
		echo "shellcheck not installed: skipping (CI runs it)"; \
	fi

test:
	$(BIN)python -m pytest --cov=artharness --cov-fail-under=$(COV_FLOOR)

fmt:
	$(BIN)ruff check --fix src tests pipeline benchmark
