# --- User Parameters


# --- Internal Parameters


# --- Targets

# Default target
all: test

test check:
	julia -e 'using Coverage; clean_folder(".");'
	julia --color=yes -e 'import Pkg; Pkg.test(;coverage=true)'
	@echo
	coverage.jl

# Maintenance
clean:
	find . -name "tmp.init-pkg.*" -exec rm -rf {} \;  # init-pkg.jl files
	julia -e 'using Coverage; clean_folder(".");'
	find . -name "*.jl.*.mem" -exec rm -f {} \;  # memory allocation tracking files

# Setup Julia
setup:
	julia --project=`pwd`/bin --startup-file=no \
		-e 'import Pkg; Pkg.instantiate()'
	julia --project=`pwd`/test --startup-file=no \
		-e 'import Pkg; Pkg.instantiate()'

# Phony Targets
.PHONY: all clean setup \
        test check
