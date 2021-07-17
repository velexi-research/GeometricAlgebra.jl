# --- User Parameters


# --- Internal Parameters


# --- Targets

# Default target
all: test

test-cmd check-cmd:
	@cd test; julia -e 'using Coverage; clean_folder("..");'
	julia --color=yes -e 'import Pkg; Pkg.test(;coverage=true)'
	coverage.jl

test check: export JULIA_TEST_FAIL_FAST = true

test check: test-cmd

test-full check-full: export JULIA_TEST_FAIL_FAST = false

test-full check-full: test-cmd

# Maintenance
clean:
	find . -name "tmp.init-pkg.*" -exec rm -rf {} \;  # init-pkg.jl files
	cd test; julia -e 'using Coverage; clean_folder("..");'
	find . -name "*.jl.*.mem" -exec rm -f {} \;  # memory allocation tracking files

# Setup Julia
setup:
	julia --project=`pwd`/bin --startup-file=no \
		-e 'import Pkg; Pkg.instantiate()'
	julia --project=`pwd`/test --startup-file=no \
		-e 'import Pkg; Pkg.instantiate()'

# Phony Targets
.PHONY: all clean setup \
        test-cmd check-cmd \
		test check \
        test-full check-full
