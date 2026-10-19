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

# ALL (default), empty (no affected tests), or test files: CI passes its
# selection through the TESTS environment variable. The recipe reads it
# quoted only ("$$TESTS"); tr + xargs split the selector's space-joined
# args format into runner arguments, so a file name is never evaluated as
# Make or shell code.
TESTS ?= ALL
export TESTS

test:
	@if [ "$$TESTS" = ALL ]; then \
		COVERAGE_CORE=$(COVERAGE_CORE) $(BIN)python -m pytest --cov=artharness --cov-fail-under=$(COV_FLOOR); \
	elif [ -z "$$TESTS" ]; then \
		echo "tests: none affected"; \
	else \
		printf '%s\n' "$$TESTS" | tr ' ' '\n' | xargs -d '\n' -r $(BIN)python scripts/ci/run_parallel.py --runner pytest; \
	fi

fmt:
	$(BIN)ruff check --fix src tests pipeline benchmark

# Whole ART workflow at toy scale; see docs/DEMO.md. Output in demo-run/WALKTHROUGH.md.
demo:
	$(BIN)python -m artharness.demo --out demo-run --force

demo-model:
	$(BIN)python -m artharness.demo --agent ollama --out demo-run --force
