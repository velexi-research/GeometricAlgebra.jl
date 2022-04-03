# --- User Parameters


# --- Internal Parameters

# Package variables
PKG_DIR=src

# Julia environment
export JULIA_PROJECT = @.

# --- Targets

# Default target
all: fast-test

# Testing
test:
	@echo Removing old coverage files
	find . -name "*.jl.*.cov" -exec rm -f {} \;
	@echo
	@echo Running tests
	julia --color=yes -e 'import Pkg; Pkg.test(; coverage=true)'
	@echo
	@echo Generating code coverage report
	@jlcoverage

fast-test: export JLTEST_FAIL_FAST=true
fast-test: test

# Code quality
codestyle:
	@echo Checking code style
	@jlcodestyle -v $(PKG_DIR)

# Documentation
docs:
	julia --project=docs --compile=min -O0 docs/make.jl

# Maintenance
clean:
	@echo Removing coverage files
	find . -name "*.jl.*.cov" -exec rm -f {} \;

spotless: clean
	find . -name "Manifest.toml" -exec rm -rf {} \;  # Manifest.toml files

# Phony Targets
.PHONY: all test fast-test codestyle docs \
		clean spotless
