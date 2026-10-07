FC      = gfortran

FFLAGS_CHECK = -O0 -g -Wall -Wextra -Wconversion-extra \
               -fcheck=all -ffpe-trap=invalid,zero,overflow \
               -fbacktrace

FFLAGS_OPT   = -O3 -march=native

BUILD_DIR = build

.PHONY: all check steady-residual friction-check docs clean

all: steady-residual

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

steady-residual: $(BUILD_DIR)
	$(FC) $(FFLAGS_CHECK) \
	    src/fortran77/steady_residual.f \
	    -o $(BUILD_DIR)/steady_residual

friction-check: $(BUILD_DIR)
	$(FC) $(FFLAGS_CHECK) \
	    src/fortran77/models/friction_swamee_jain.f \
	    tests/fortran77/friction_check.f \
	    -o $(BUILD_DIR)/friction_check

check: steady-residual friction-check
	./$(BUILD_DIR)/steady_residual
	./$(BUILD_DIR)/friction_check

docs:
	cd docs && pdflatex -halt-on-error model.tex
	cd docs && pdflatex -halt-on-error model.tex
	cd docs && pdflatex -halt-on-error model.tex

clean:
	rm -rf $(BUILD_DIR)
	rm -f docs/model.aux docs/model.log docs/model.out \
	      docs/model.toc docs/model.fls docs/model.fdb_latexmk \
	      docs/model.synctex.gz
