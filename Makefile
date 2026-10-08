FC      = gfortran

FFLAGS_CHECK = -O0 -g -Wall -Wextra -Wconversion-extra \
               -fcheck=all -ffpe-trap=invalid,zero,overflow \
               -fbacktrace

FFLAGS_OPT   = -O3 -march=native

BUILD_DIR = build

.PHONY: m9d2-temporal-check all check steady-residual friction-check jacobian-check newton-check timestep-check integration-check adaptive-check adaptive-schedule-check adaptive-schedule-pressure-check input-validation-check characterize dynamics waves wave-analysis modal-analysis small-verify docs clean m9e2-interface-check

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
	    src/fortran77/input_validation.f \
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
	    src/fortran77/input_validation.f \
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
	    src/fortran77/input_validation.f \
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
	    src/fortran77/input_validation.f \
	    src/fortran77/newton_solver.f \
	    src/fortran77/transient_step.f \
	    src/fortran77/transient_integrate.f \
	    src/fortran77/transient_adaptive.f \
	    tests/fortran77/adaptive_integration_check.f \
	    -o $(BUILD_DIR)/adaptive_integration_check

adaptive-schedule-check adaptive-schedule-pressure-check: $(BUILD_DIR)
	$(FC) $(FFLAGS_CHECK) \
	    src/fortran77/models/eos_constant_z.f \
	    src/fortran77/models/friction_swamee_jain.f \
	    src/fortran77/transient_residual.f \
	    src/fortran77/transient_jacobian.f \
	    src/fortran77/linear_solve.f \
	    src/fortran77/input_validation.f \
	    src/fortran77/newton_solver.f \
	    src/fortran77/transient_step.f \
	    src/fortran77/transient_integrate.f \
	    src/fortran77/transient_adaptive_schedule.f \
	    tests/fortran77/$(subst -,_,$@).f \
	    -o $(BUILD_DIR)/$@

m9b3-negative-check: $(BUILD_DIR)
	$(FC) $(FFLAGS_CHECK) \
	    src/fortran77/models/eos_constant_z.f \
	    src/fortran77/models/friction_swamee_jain.f \
	    src/fortran77/transient_residual.f \
	    src/fortran77/transient_jacobian.f \
	    src/fortran77/linear_solve.f \
	    src/fortran77/input_validation.f \
	    src/fortran77/newton_solver.f \
	    src/fortran77/transient_step.f \
	    src/fortran77/transient_integrate.f \
	    src/fortran77/transient_adaptive.f \
	    src/fortran77/transient_adaptive_schedule.f \
	    tests/fortran77/m9b3_negative_check.f \
	    -o $(BUILD_DIR)/m9b3_negative_check

input-validation-check: $(BUILD_DIR)
	$(FC) $(FFLAGS_CHECK) \
	    src/fortran77/models/eos_constant_z.f \
	    src/fortran77/models/friction_swamee_jain.f \
	    src/fortran77/transient_residual.f \
	    src/fortran77/transient_jacobian.f \
	    src/fortran77/linear_solve.f \
	    src/fortran77/input_validation.f \
	    src/fortran77/newton_solver.f \
	    src/fortran77/transient_step.f \
	    tests/fortran77/input_validation_check.f \
	    -o $(BUILD_DIR)/input_validation_check


m9c2a-failure-check: $(BUILD_DIR)
	$(FC) $(FFLAGS_CHECK) \
	    src/fortran77/models/eos_constant_z.f \
	    src/fortran77/models/friction_swamee_jain.f \
	    src/fortran77/transient_residual.f \
	    src/fortran77/transient_jacobian.f \
	    src/fortran77/linear_solve.f \
	    src/fortran77/input_validation.f \
	    src/fortran77/newton_solver.f \
	    src/fortran77/transient_step.f \
	    src/fortran77/transient_integrate.f \
	    src/fortran77/transient_adaptive.f \
	    src/fortran77/transient_adaptive_schedule.f \
	    tests/fortran77/m9c2a_failure_contract_check.f \
	    -o $(BUILD_DIR)/m9c2a_failure_contract_check


m9c2b-failure-check: $(BUILD_DIR)
	$(FC) $(FFLAGS_CHECK) \
	    src/fortran77/models/eos_constant_z.f \
	    src/fortran77/models/friction_swamee_jain.f \
	    src/fortran77/transient_residual.f \
	    src/fortran77/transient_jacobian.f \
	    src/fortran77/linear_solve.f \
	    src/fortran77/input_validation.f \
	    src/fortran77/newton_solver.f \
	    src/fortran77/transient_step.f \
	    src/fortran77/transient_integrate.f \
	    src/fortran77/transient_adaptive.f \
	    src/fortran77/transient_adaptive_schedule.f \
	    tests/fortran77/m9c2b_failure_contract_check.f \
	    -o $(BUILD_DIR)/m9c2b_failure_contract_check

m9d2-temporal-check: $(BUILD_DIR)
	$(FC) $(FFLAGS_CHECK) \
	    src/fortran77/models/eos_constant_z.f \
	    src/fortran77/models/friction_swamee_jain.f \
	    src/fortran77/transient_residual.f \
	    src/fortran77/transient_jacobian.f \
	    src/fortran77/linear_solve.f \
	    src/fortran77/input_validation.f \
	    src/fortran77/newton_solver.f \
	    src/fortran77/transient_step.f \
	    src/fortran77/transient_integrate.f \
	    tests/fortran77/m9d2_temporal.f90 \
	    -o $(BUILD_DIR)/m9d2_temporal

m9d3-spatial-check: $(BUILD_DIR)
	$(FC) $(FFLAGS_CHECK) \
	    src/fortran77/models/eos_constant_z.f \
	    src/fortran77/models/friction_swamee_jain.f \
	    src/fortran77/transient_residual.f \
	    src/fortran77/transient_jacobian.f \
	    src/fortran77/linear_solve.f \
	    src/fortran77/input_validation.f \
	    src/fortran77/newton_solver.f \
	    src/fortran77/transient_step.f \
	    src/fortran77/transient_integrate.f \
	    tests/fortran77/m9d3_spatial_regression.f90 \
	    -o $(BUILD_DIR)/m9d3_spatial_regression

m9d4-inventory-check: $(BUILD_DIR)
	$(FC) $(FFLAGS_CHECK) \
	    src/fortran77/models/eos_constant_z.f \
	    src/fortran77/models/friction_swamee_jain.f \
	    src/fortran77/transient_residual.f \
	    src/fortran77/transient_jacobian.f \
	    src/fortran77/linear_solve.f \
	    src/fortran77/input_validation.f \
	    src/fortran77/newton_solver.f \
	    src/fortran77/transient_step.f \
	    src/fortran77/transient_integrate.f \
	    tests/fortran77/m9d4_inventory_regression.f90 \
	    -o $(BUILD_DIR)/m9d4_inventory_regression

m9d5-halfcell-check: $(BUILD_DIR)
	$(FC) $(FFLAGS_CHECK) \
	    src/fortran77/models/eos_constant_z.f \
	    src/fortran77/models/friction_swamee_jain.f \
	    src/fortran77/transient_residual.f \
	    src/fortran77/transient_jacobian.f \
	    tests/fortran77/m9d5_halfcell_regression.f90 \
	    -o $(BUILD_DIR)/m9d5_halfcell_regression

m9e2-interface-check: $(BUILD_DIR)
	mkdir -p $(BUILD_DIR)/interfaces
	$(FC) $(FFLAGS_CHECK) -J $(BUILD_DIR)/interfaces \
	    -c src/interfaces/pipe_solver_api.f90 \
	    -o $(BUILD_DIR)/interfaces/pipe_solver_api.o
	$(FC) $(FFLAGS_CHECK) -I $(BUILD_DIR)/interfaces \
	    $(BUILD_DIR)/interfaces/pipe_solver_api.o \
	    src/fortran77/models/eos_constant_z.f \
	    src/fortran77/models/friction_swamee_jain.f \
	    src/fortran77/transient_residual.f \
	    src/fortran77/transient_jacobian.f \
	    src/fortran77/linear_solve.f \
	    src/fortran77/input_validation.f \
	    src/fortran77/newton_solver.f \
	    src/fortran77/transient_step.f \
	    tests/fortran77/m9e2_interface_check.f90 \
	    -o $(BUILD_DIR)/m9e2_interface_check

check: steady-residual friction-check jacobian-check newton-check timestep-check integration-check adaptive-check adaptive-schedule-check adaptive-schedule-pressure-check input-validation-check m9b3-negative-check m9c2a-failure-check m9c2b-failure-check m9c3-partial-failure-check m9d2-temporal-check m9d3-spatial-check m9d4-inventory-check m9d5-halfcell-check m9e2-interface-check
	./$(BUILD_DIR)/steady_residual
	./$(BUILD_DIR)/friction_check
	./$(BUILD_DIR)/jacobian_check
	./$(BUILD_DIR)/newton_check
	./$(BUILD_DIR)/timestep_check
	./$(BUILD_DIR)/transient_integration_check
	./$(BUILD_DIR)/adaptive_integration_check
	./$(BUILD_DIR)/adaptive-schedule-check
	./$(BUILD_DIR)/adaptive-schedule-pressure-check
	./$(BUILD_DIR)/input_validation_check
	./$(BUILD_DIR)/m9b3_negative_check
	./$(BUILD_DIR)/m9c2a_failure_contract_check
	./$(BUILD_DIR)/m9c2b_failure_contract_check
	./$(BUILD_DIR)/m9c3_partial_failure_check
	./$(BUILD_DIR)/m9d2_temporal
	./$(BUILD_DIR)/m9d3_spatial_regression
	./$(BUILD_DIR)/m9d4_inventory_regression
	./$(BUILD_DIR)/m9d5_halfcell_regression
	./$(BUILD_DIR)/m9e2_interface_check

characterize: $(BUILD_DIR)
	$(FC) $(FFLAGS_OPT) -Wall -Wextra -Wconversion-extra \
	    src/fortran77/models/eos_constant_z.f \
	    src/fortran77/models/friction_swamee_jain.f \
	    src/fortran77/transient_residual.f \
	    src/fortran77/transient_jacobian.f \
	    src/fortran77/linear_solve.f \
	    src/fortran77/input_validation.f \
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
	    src/fortran77/input_validation.f \
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
	    src/fortran77/input_validation.f \
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
	    src/fortran77/input_validation.f \
	    src/fortran77/newton_solver.f \
	    src/fortran77/transient_step.f \
	    src/fortran77/transient_integrate.f \
	    tests/fortran77/small_perturbation.f \
	    -o $(BUILD_DIR)/small_perturbation
	./$(BUILD_DIR)/small_perturbation
	python3 benchmarks/verify_small_m7b.py

m9c3-partial-failure-check: $(BUILD_DIR)
	$(FC) $(FFLAGS_CHECK) -Wno-unused-dummy-argument \
	    src/fortran77/models/eos_constant_z.f \
	    src/fortran77/models/friction_swamee_jain.f \
	    src/fortran77/transient_residual.f \
	    src/fortran77/input_validation.f \
	    src/fortran77/transient_integrate.f \
	    src/fortran77/transient_adaptive.f \
	    src/fortran77/transient_adaptive_schedule.f \
	    tests/fortran77/m9c3_partial_failure_check.f \
	    -o $(BUILD_DIR)/m9c3_partial_failure_check
