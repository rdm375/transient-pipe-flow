FC      = gfortran

FFLAGS_CHECK = -O0 -g -Wall -Wextra -Wconversion-extra \
               -fcheck=all -ffpe-trap=invalid,zero,overflow \
               -fbacktrace

FFLAGS_OPT   = -O3 -march=native

BUILD_DIR = build

.PHONY: all check steady-residual friction-check jacobian-check newton-check timestep-check docs clean

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


jacobian-check: $(BUILD_DIR)
	$(FC) $(FFLAGS_CHECK) \
	    src/fortran77/models/eos_constant_z.f \
	    src/fortran77/models/friction_swamee_jain.f \
	    src/fortran77/transient_residual.f \
	    src/fortran77/transient_jacobian.f \
	    tests/fortran77/jacobian_check.f \
	    -o $(BUILD_DIR)/jacobian_check


newton-check: $(BUILD_DIR)
	$(FC) $(FFLAGS_CHECK) \
	    src/fortran77/models/eos_constant_z.f \
	    src/fortran77/models/friction_swamee_jain.f \
	    src/fortran77/transient_residual.f \
	    src/fortran77/transient_jacobian.f \
	    src/fortran77/linear_solve.f \
	    src/fortran77/newton_solver.f \
	    tests/fortran77/newton_check.f \
	    -o $(BUILD_DIR)/newton_check


timestep-check: $(BUILD_DIR)
	$(FC) $(FFLAGS_CHECK) \
	    src/fortran77/models/eos_constant_z.f \
	    src/fortran77/models/friction_swamee_jain.f \
	    src/fortran77/transient_residual.f \
	    src/fortran77/transient_jacobian.f \
	    src/fortran77/linear_solve.f \
	    src/fortran77/newton_solver.f \
	    src/fortran77/transient_step.f \
	    tests/fortran77/timestep_check.f \
	    -o $(BUILD_DIR)/timestep_check

check: steady-residual friction-check jacobian-check newton-check timestep-check
	./$(BUILD_DIR)/steady_residual
	./$(BUILD_DIR)/friction_check
	./$(BUILD_DIR)/jacobian_check
	./$(BUILD_DIR)/newton_check
	./$(BUILD_DIR)/timestep_check

docs:
	cd docs && pdflatex -halt-on-error model.tex
	cd docs && pdflatex -halt-on-error model.tex
	cd docs && pdflatex -halt-on-error model.tex

clean:
	rm -rf $(BUILD_DIR)
	rm -f docs/model.aux docs/model.log docs/model.out \
	      docs/model.toc docs/model.fls docs/model.fdb_latexmk \
	      docs/model.synctex.gz
