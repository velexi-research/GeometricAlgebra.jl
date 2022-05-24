# --- Makefile parameters

# Default rule
.DEFAULT_GOAL := fast-test

# Package variables
PKG_DIR=src

# Julia environment
export JULIA_PROJECT = @.

# --- Testing rules

.PHONY: test fast-test

## Run all tests
test:
	@echo Removing old coverage files
	find . -name "*.jl.*.cov" -exec rm -f {} \;
	@echo
	@echo Running tests
	julia --color=yes -e 'import Pkg; Pkg.test(; coverage=true)'
	@echo
	@echo Generating code coverage report
	@jlcoverage

## Run tests in fail-fast mode (i.e., stop at first failure)
fast-test: export JLTEST_FAIL_FAST=true
fast-test: test

# --- Code quality rules

.PHONY: codestyle

## Check codestyle
codestyle:
	@echo Checking code style
	@jlcodestyle -v $(PKG_DIR)

# --- Documentation rules

.PHONY: docs

## Generate package documentation.
docs:
	julia --project=docs --compile=min -O0 docs/make.jl

# --- Utility rules

.PHONY: clean

## Remove files and directories automatically generated during development (e.g., coverage
## files).
clean:
	@echo Removing coverage files
	find . -name "*.jl.*.cov" -exec rm -f {} \;

## Remove files and directories automatically generated during development (e.g., coverage
## files) and project setup (e.g., `Manifest.toml` files).
spotless: clean
	@echo Removing Manifest.toml files
	find . -name "Manifest.toml" -exec rm -rf {} \;

# --- Makefile Self-Documentation

# Inspired by
# <http://marmelab.com/blog/2016/02/29/auto-documented-makefile.html>
#
# sed script explained:
# /^##/:
# 	* save line in hold space
# 	* purge line
# 	* Loop:
# 		* append newline + line to hold space
# 		* go to next line
# 		* if line starts with doc comment, strip comment character off and loop
# 	* remove target prerequisites
# 	* append hold space (+ newline) to line
# 	* replace newline plus comments by `---`
# 	* print line
# Separate expressions are necessary because labels cannot be delimited by
# semicolon; see <http://stackoverflow.com/a/11799865/1968>

.PHONY: help

## Display this list of available rules
help:
	@echo "$$(tput bold)Default rule:$$(tput sgr0) ${.DEFAULT_GOAL}"
	@echo
	@echo "$$(tput bold)Available rules:$$(tput sgr0)"
	@echo
	@sed -n -e "/^## / { \
		h; \
		s/.*//; \
		:doc" \
		-e "H; \
		n; \
		s/^## //; \
		t doc" \
		-e "s/:.*//; \
		G; \
		s/\\n## /---/; \
		s/\\n/ /g; \
		p; \
	}" ${MAKEFILE_LIST} \
	| LC_ALL='C' sort --ignore-case \
	| awk -F '---' \
		-v ncol=$$(tput cols) \
		-v indent=19 \
		-v col_on="$$(tput setaf 6)" \
		-v col_off="$$(tput sgr0)" \
	'{ \
		printf "%s%*s%s ", col_on, -indent, $$1, col_off; \
		n = split($$2, words, " "); \
		line_length = ncol - indent; \
		for (i = 1; i <= n; i++) { \
			line_length -= length(words[i]) + 1; \
			if (line_length <= 0) { \
				line_length = ncol - indent - length(words[i]) - 1; \
				printf "\n%*s ", -indent, " "; \
			} \
			printf "%s ", words[i]; \
		} \
		printf "\n"; \
	}' \
	| more $(shell test $(shell uname) = Darwin && echo '--no-init --raw-control-chars')
