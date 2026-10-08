FC      = gfortran

FFLAGS_CHECK = -O0 -g -Wall -Wextra -Wconversion-extra \
               -fcheck=all -ffpe-trap=invalid,zero,overflow \
               -fbacktrace

FFLAGS_OPT   = -O3 -march=native

BUILD_DIR = build

.PHONY: all check steady-residual friction-check jacobian-check newton-check timestep-check integration-check adaptive-check characterize dynamics waves wave-analysis modal-analysis small-verify docs clean

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


integration-check: $(BUILD_DIR)
	$(FC) $(FFLAGS_CHECK) \
	    src/fortran77/models/eos_constant_z.f \
	    src/fortran77/models/friction_swamee_jain.f \
	    src/fortran77/transient_residual.f \
	    src/fortran77/transient_jacobian.f \
	    src/fortran77/linear_solve.f \
	    src/fortran77/newton_solver.f \
	    src/fortran77/transient_step.f \
	    src/fortran77/transient_integrate.f \
	    tests/fortran77/transient_integration_check.f \
	    -o $(BUILD_DIR)/transient_integration_check

adaptive-check: $(BUILD_DIR)
	$(FC) $(FFLAGS_CHECK) \
	    src/fortran77/models/eos_constant_z.f \
	    src/fortran77/models/friction_swamee_jain.f \
	    src/fortran77/transient_residual.f \
	    src/fortran77/transient_jacobian.f \
	    src/fortran77/linear_solve.f \
	    src/fortran77/newton_solver.f \
	    src/fortran77/transient_step.f \
	    src/fortran77/transient_integrate.f \
	    src/fortran77/transient_adaptive.f \
	    tests/fortran77/adaptive_integration_check.f \
	    -o $(BUILD_DIR)/adaptive_integration_check

check: steady-residual friction-check jacobian-check newton-check timestep-check integration-check adaptive-check
	./$(BUILD_DIR)/steady_residual
	./$(BUILD_DIR)/friction_check
	./$(BUILD_DIR)/jacobian_check
	./$(BUILD_DIR)/newton_check
	./$(BUILD_DIR)/timestep_check
	./$(BUILD_DIR)/transient_integration_check
	./$(BUILD_DIR)/adaptive_integration_check

characterize: $(BUILD_DIR)
	$(FC) $(FFLAGS_OPT) -Wall -Wextra -Wconversion-extra \
	    src/fortran77/models/eos_constant_z.f \
	    src/fortran77/models/friction_swamee_jain.f \
	    src/fortran77/transient_residual.f \
	    src/fortran77/transient_jacobian.f \
	    src/fortran77/linear_solve.f \
	    src/fortran77/newton_solver.f \
	    src/fortran77/transient_step.f \
	    src/fortran77/transient_integrate.f \
	    tests/fortran77/numerical_characterization.f \
	    -o $(BUILD_DIR)/numerical_characterization
	./$(BUILD_DIR)/numerical_characterization

dynamics: $(BUILD_DIR)
	$(FC) $(FFLAGS_OPT) -Wall -Wextra -Wconversion-extra \
	    src/fortran77/models/eos_constant_z.f \
	    src/fortran77/models/friction_swamee_jain.f \
	    src/fortran77/transient_residual.f \
	    src/fortran77/transient_jacobian.f \
	    src/fortran77/linear_solve.f \
	    src/fortran77/newton_solver.f \
	    src/fortran77/transient_step.f \
	    src/fortran77/transient_integrate.f \
	    tests/fortran77/dynamic_characterization.f \
	    -o $(BUILD_DIR)/dynamic_characterization
	./$(BUILD_DIR)/dynamic_characterization

waves: $(BUILD_DIR)
	$(FC) $(FFLAGS_OPT) -Wall -Wextra -Wconversion-extra \
	    src/fortran77/models/eos_constant_z.f \
	    src/fortran77/models/friction_swamee_jain.f \
	    src/fortran77/transient_residual.f \
	    src/fortran77/transient_jacobian.f \
	    src/fortran77/linear_solve.f \
	    src/fortran77/newton_solver.f \
	    src/fortran77/transient_step.f \
	    src/fortran77/transient_integrate.f \
	    tests/fortran77/wave_characterization.f \
	    -o $(BUILD_DIR)/wave_characterization
	./$(BUILD_DIR)/wave_characterization

wave-analysis: waves
	python3 benchmarks/analyze_m7b.py

docs:
	cd docs && pdflatex -halt-on-error model.tex
	cd docs && pdflatex -halt-on-error model.tex
	cd docs && pdflatex -halt-on-error model.tex

clean:
	rm -rf $(BUILD_DIR)
	rm -f docs/model.aux docs/model.log docs/model.out \
	      docs/model.toc docs/model.fls docs/model.fdb_latexmk \
	      docs/model.synctex.gz

modal-analysis: wave-analysis
	$(FC) $(FFLAGS_OPT) -Wall -Wextra -Wconversion-extra \
	    src/fortran77/models/eos_constant_z.f \
	    src/fortran77/models/friction_swamee_jain.f \
	    src/fortran77/transient_jacobian.f \
	    tests/fortran77/modal_jacobian.f \
	    -o $(BUILD_DIR)/modal_jacobian
	./$(BUILD_DIR)/modal_jacobian
	python3 benchmarks/modal_m7b.py

# Small-perturbation verification uses the production nonlinear stepper.
small-verify: modal-analysis
	$(FC) $(FFLAGS_OPT) -Wall -Wextra -Wconversion-extra \
	    src/fortran77/models/eos_constant_z.f \
	    src/fortran77/models/friction_swamee_jain.f \
	    src/fortran77/transient_residual.f \
	    src/fortran77/transient_jacobian.f \
	    src/fortran77/linear_solve.f \
	    src/fortran77/newton_solver.f \
	    src/fortran77/transient_step.f \
	    src/fortran77/transient_integrate.f \
	    tests/fortran77/small_perturbation.f \
	    -o $(BUILD_DIR)/small_perturbation
	./$(BUILD_DIR)/small_perturbation
	python3 benchmarks/verify_small_m7b.py
