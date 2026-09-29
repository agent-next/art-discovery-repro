.PHONY: setup check lint shellcheck test fmt demo demo-model

# Use the project venv when it exists (make setup / uv sync); fall back to PATH (CI).
BIN := $(if $(wildcard .venv/bin/ruff),.venv/bin/,)
COV_FLOOR := 94
# sys.monitoring (3.12+) traces the numeric demo/array code ~2.5x faster than ctrace.
COVERAGE_CORE := $(shell $(BIN)python -c 'import sys; print("sysmon" if sys.version_info >= (3, 12) else "ctrace")')

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
	COVERAGE_CORE=$(COVERAGE_CORE) $(BIN)python -m pytest --cov=artharness --cov-fail-under=$(COV_FLOOR)

fmt:
	$(BIN)ruff check --fix src tests pipeline benchmark

# Whole ART workflow at toy scale; see docs/DEMO.md. Output in demo-run/WALKTHROUGH.md.
demo:
	$(BIN)python -m artharness.demo --out demo-run --force

demo-model:
	$(BIN)python -m artharness.demo --agent ollama --out demo-run --force
