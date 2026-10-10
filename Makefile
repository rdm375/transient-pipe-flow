build/obj/check/tests/fortran77/m9c3_partial_failure_check.o: private FFLAGS_CHECK += -Wno-unused-dummy-argument
FC      = gfortran

FFLAGS_CHECK = -O0 -g -Wall -Wextra -Wconversion-extra \
               -fcheck=all -ffpe-trap=invalid,zero,overflow \
               -fbacktrace

FFLAGS_OPT   = -O3 -march=native

BUILD_DIR = build

.PHONY: m9e3-failure-check m9d2-temporal-check all check steady-residual friction-check jacobian-check newton-check timestep-check integration-check adaptive-check adaptive-schedule-check adaptive-schedule-pressure-check input-validation-check characterize dynamics waves wave-analysis modal-analysis small-verify docs clean m9e2-interface-check m9e2b-integration-check m9e2c-adaptive-check m9e2d-schedule-check

all: steady-residual

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

steady-residual: $(BUILD_DIR)/steady_residual

$(BUILD_DIR)/steady_residual: $(BUILD_DIR)/obj/check/src/fortran77/steady_residual.o
	@mkdir -p $(@D)
	$(FC) $(FFLAGS_CHECK) -I $(BUILD_DIR)/interfaces $^ -o $@

friction-check: $(BUILD_DIR)/friction_check

$(BUILD_DIR)/friction_check: $(BUILD_DIR)/obj/check/src/fortran77/models/friction_swamee_jain.o $(BUILD_DIR)/obj/check/tests/fortran77/friction_check.o
	@mkdir -p $(@D)
	$(FC) $(FFLAGS_CHECK) -I $(BUILD_DIR)/interfaces $^ -o $@


jacobian-check: $(BUILD_DIR)/jacobian_check

$(BUILD_DIR)/jacobian_check: $(BUILD_DIR)/obj/check/src/fortran77/models/eos_constant_z.o $(BUILD_DIR)/obj/check/src/fortran77/models/friction_swamee_jain.o $(BUILD_DIR)/obj/check/src/fortran77/transient_residual.o $(BUILD_DIR)/obj/check/src/fortran77/transient_jacobian.o $(BUILD_DIR)/obj/check/tests/fortran77/jacobian_check.o
	@mkdir -p $(@D)
	$(FC) $(FFLAGS_CHECK) -I $(BUILD_DIR)/interfaces $^ -o $@


newton-check: $(BUILD_DIR)/newton_check

$(BUILD_DIR)/newton_check: $(BUILD_DIR)/obj/check/src/fortran77/models/eos_constant_z.o $(BUILD_DIR)/obj/check/src/fortran77/models/friction_swamee_jain.o $(BUILD_DIR)/obj/check/src/fortran77/transient_residual.o $(BUILD_DIR)/obj/check/src/fortran77/transient_jacobian.o $(BUILD_DIR)/obj/check/src/fortran77/linear_solve.o $(BUILD_DIR)/obj/check/src/fortran77/input_validation.o $(BUILD_DIR)/obj/check/src/fortran77/newton_solver.o $(BUILD_DIR)/obj/check/tests/fortran77/newton_check.o
	@mkdir -p $(@D)
	$(FC) $(FFLAGS_CHECK) -I $(BUILD_DIR)/interfaces $^ -o $@


timestep-check: $(BUILD_DIR)/timestep_check

$(BUILD_DIR)/timestep_check: $(BUILD_DIR)/obj/check/src/fortran77/models/eos_constant_z.o $(BUILD_DIR)/obj/check/src/fortran77/models/friction_swamee_jain.o $(BUILD_DIR)/obj/check/src/fortran77/transient_residual.o $(BUILD_DIR)/obj/check/src/fortran77/transient_jacobian.o $(BUILD_DIR)/obj/check/src/fortran77/linear_solve.o $(BUILD_DIR)/obj/check/src/fortran77/input_validation.o $(BUILD_DIR)/obj/check/src/fortran77/newton_solver.o $(BUILD_DIR)/obj/check/src/fortran77/transient_step.o $(BUILD_DIR)/obj/check/tests/fortran77/timestep_check.o
	@mkdir -p $(@D)
	$(FC) $(FFLAGS_CHECK) -I $(BUILD_DIR)/interfaces $^ -o $@


integration-check: $(BUILD_DIR)/transient_integration_check

$(BUILD_DIR)/transient_integration_check: $(BUILD_DIR)/obj/check/src/fortran77/models/eos_constant_z.o $(BUILD_DIR)/obj/check/src/fortran77/models/friction_swamee_jain.o $(BUILD_DIR)/obj/check/src/fortran77/transient_residual.o $(BUILD_DIR)/obj/check/src/fortran77/transient_jacobian.o $(BUILD_DIR)/obj/check/src/fortran77/linear_solve.o $(BUILD_DIR)/obj/check/src/fortran77/input_validation.o $(BUILD_DIR)/obj/check/src/fortran77/newton_solver.o $(BUILD_DIR)/obj/check/src/fortran77/transient_step.o $(BUILD_DIR)/obj/check/src/fortran77/transient_integrate.o $(BUILD_DIR)/obj/check/tests/fortran77/transient_integration_check.o
	@mkdir -p $(@D)
	$(FC) $(FFLAGS_CHECK) -I $(BUILD_DIR)/interfaces $^ -o $@

adaptive-check: $(BUILD_DIR)/adaptive_integration_check

$(BUILD_DIR)/adaptive_integration_check: $(BUILD_DIR)/obj/check/src/fortran77/models/eos_constant_z.o $(BUILD_DIR)/obj/check/src/fortran77/models/friction_swamee_jain.o $(BUILD_DIR)/obj/check/src/fortran77/transient_residual.o $(BUILD_DIR)/obj/check/src/fortran77/transient_jacobian.o $(BUILD_DIR)/obj/check/src/fortran77/linear_solve.o $(BUILD_DIR)/obj/check/src/fortran77/input_validation.o $(BUILD_DIR)/obj/check/src/fortran77/newton_solver.o $(BUILD_DIR)/obj/check/src/fortran77/transient_step.o $(BUILD_DIR)/obj/check/src/fortran77/transient_integrate.o $(BUILD_DIR)/obj/check/src/fortran77/transient_adaptive.o $(BUILD_DIR)/obj/check/tests/fortran77/adaptive_integration_check.o
	@mkdir -p $(@D)
	$(FC) $(FFLAGS_CHECK) -I $(BUILD_DIR)/interfaces $^ -o $@

adaptive-schedule-check: $(BUILD_DIR)/adaptive-schedule-check

$(BUILD_DIR)/adaptive-schedule-check: $(BUILD_DIR)/obj/check/src/fortran77/models/eos_constant_z.o $(BUILD_DIR)/obj/check/src/fortran77/models/friction_swamee_jain.o $(BUILD_DIR)/obj/check/src/fortran77/transient_residual.o $(BUILD_DIR)/obj/check/src/fortran77/transient_jacobian.o $(BUILD_DIR)/obj/check/src/fortran77/linear_solve.o $(BUILD_DIR)/obj/check/src/fortran77/input_validation.o $(BUILD_DIR)/obj/check/src/fortran77/newton_solver.o $(BUILD_DIR)/obj/check/src/fortran77/transient_step.o $(BUILD_DIR)/obj/check/src/fortran77/transient_integrate.o $(BUILD_DIR)/obj/check/src/fortran77/transient_adaptive_schedule.o $(BUILD_DIR)/obj/check/tests/fortran77/adaptive_schedule_check.o
	$(FC) $(FFLAGS_CHECK) $^ -o $@

adaptive-schedule-pressure-check: $(BUILD_DIR)/adaptive-schedule-pressure-check

$(BUILD_DIR)/adaptive-schedule-pressure-check: $(BUILD_DIR)/obj/check/src/fortran77/models/eos_constant_z.o $(BUILD_DIR)/obj/check/src/fortran77/models/friction_swamee_jain.o $(BUILD_DIR)/obj/check/src/fortran77/transient_residual.o $(BUILD_DIR)/obj/check/src/fortran77/transient_jacobian.o $(BUILD_DIR)/obj/check/src/fortran77/linear_solve.o $(BUILD_DIR)/obj/check/src/fortran77/input_validation.o $(BUILD_DIR)/obj/check/src/fortran77/newton_solver.o $(BUILD_DIR)/obj/check/src/fortran77/transient_step.o $(BUILD_DIR)/obj/check/src/fortran77/transient_integrate.o $(BUILD_DIR)/obj/check/src/fortran77/transient_adaptive_schedule.o $(BUILD_DIR)/obj/check/tests/fortran77/adaptive_schedule_pressure_check.o
	$(FC) $(FFLAGS_CHECK) $^ -o $@


m9b3-negative-check: $(BUILD_DIR)/m9b3_negative_check

$(BUILD_DIR)/m9b3_negative_check: $(BUILD_DIR)/obj/check/src/fortran77/models/eos_constant_z.o $(BUILD_DIR)/obj/check/src/fortran77/models/friction_swamee_jain.o $(BUILD_DIR)/obj/check/src/fortran77/transient_residual.o $(BUILD_DIR)/obj/check/src/fortran77/transient_jacobian.o $(BUILD_DIR)/obj/check/src/fortran77/linear_solve.o $(BUILD_DIR)/obj/check/src/fortran77/input_validation.o $(BUILD_DIR)/obj/check/src/fortran77/newton_solver.o $(BUILD_DIR)/obj/check/src/fortran77/transient_step.o $(BUILD_DIR)/obj/check/src/fortran77/transient_integrate.o $(BUILD_DIR)/obj/check/src/fortran77/transient_adaptive.o $(BUILD_DIR)/obj/check/src/fortran77/transient_adaptive_schedule.o $(BUILD_DIR)/obj/check/tests/fortran77/m9b3_negative_check.o
	@mkdir -p $(@D)
	$(FC) $(FFLAGS_CHECK) -I $(BUILD_DIR)/interfaces $^ -o $@

input-validation-check: $(BUILD_DIR)/input_validation_check

$(BUILD_DIR)/input_validation_check: $(BUILD_DIR)/obj/check/src/fortran77/models/eos_constant_z.o $(BUILD_DIR)/obj/check/src/fortran77/models/friction_swamee_jain.o $(BUILD_DIR)/obj/check/src/fortran77/transient_residual.o $(BUILD_DIR)/obj/check/src/fortran77/transient_jacobian.o $(BUILD_DIR)/obj/check/src/fortran77/linear_solve.o $(BUILD_DIR)/obj/check/src/fortran77/input_validation.o $(BUILD_DIR)/obj/check/src/fortran77/newton_solver.o $(BUILD_DIR)/obj/check/src/fortran77/transient_step.o $(BUILD_DIR)/obj/check/tests/fortran77/input_validation_check.o
	@mkdir -p $(@D)
	$(FC) $(FFLAGS_CHECK) -I $(BUILD_DIR)/interfaces $^ -o $@


m9c2a-failure-check: $(BUILD_DIR)/m9c2a_failure_contract_check

$(BUILD_DIR)/m9c2a_failure_contract_check: $(BUILD_DIR)/obj/check/src/fortran77/models/eos_constant_z.o $(BUILD_DIR)/obj/check/src/fortran77/models/friction_swamee_jain.o $(BUILD_DIR)/obj/check/src/fortran77/transient_residual.o $(BUILD_DIR)/obj/check/src/fortran77/transient_jacobian.o $(BUILD_DIR)/obj/check/src/fortran77/linear_solve.o $(BUILD_DIR)/obj/check/src/fortran77/input_validation.o $(BUILD_DIR)/obj/check/src/fortran77/newton_solver.o $(BUILD_DIR)/obj/check/src/fortran77/transient_step.o $(BUILD_DIR)/obj/check/src/fortran77/transient_integrate.o $(BUILD_DIR)/obj/check/src/fortran77/transient_adaptive.o $(BUILD_DIR)/obj/check/src/fortran77/transient_adaptive_schedule.o $(BUILD_DIR)/obj/check/tests/fortran77/m9c2a_failure_contract_check.o
	@mkdir -p $(@D)
	$(FC) $(FFLAGS_CHECK) -I $(BUILD_DIR)/interfaces $^ -o $@


m9c2b-failure-check: $(BUILD_DIR)/m9c2b_failure_contract_check

$(BUILD_DIR)/m9c2b_failure_contract_check: $(BUILD_DIR)/obj/check/src/fortran77/models/eos_constant_z.o $(BUILD_DIR)/obj/check/src/fortran77/models/friction_swamee_jain.o $(BUILD_DIR)/obj/check/src/fortran77/transient_residual.o $(BUILD_DIR)/obj/check/src/fortran77/transient_jacobian.o $(BUILD_DIR)/obj/check/src/fortran77/linear_solve.o $(BUILD_DIR)/obj/check/src/fortran77/input_validation.o $(BUILD_DIR)/obj/check/src/fortran77/newton_solver.o $(BUILD_DIR)/obj/check/src/fortran77/transient_step.o $(BUILD_DIR)/obj/check/src/fortran77/transient_integrate.o $(BUILD_DIR)/obj/check/src/fortran77/transient_adaptive.o $(BUILD_DIR)/obj/check/src/fortran77/transient_adaptive_schedule.o $(BUILD_DIR)/obj/check/tests/fortran77/m9c2b_failure_contract_check.o
	@mkdir -p $(@D)
	$(FC) $(FFLAGS_CHECK) -I $(BUILD_DIR)/interfaces $^ -o $@

m9d2-temporal-check: $(BUILD_DIR)/m9d2_temporal

$(BUILD_DIR)/m9d2_temporal: $(BUILD_DIR)/obj/check/src/fortran77/models/eos_constant_z.o $(BUILD_DIR)/obj/check/src/fortran77/models/friction_swamee_jain.o $(BUILD_DIR)/obj/check/src/fortran77/transient_residual.o $(BUILD_DIR)/obj/check/src/fortran77/transient_jacobian.o $(BUILD_DIR)/obj/check/src/fortran77/linear_solve.o $(BUILD_DIR)/obj/check/src/fortran77/input_validation.o $(BUILD_DIR)/obj/check/src/fortran77/newton_solver.o $(BUILD_DIR)/obj/check/src/fortran77/transient_step.o $(BUILD_DIR)/obj/check/src/fortran77/transient_integrate.o $(BUILD_DIR)/obj/check/tests/fortran77/m9d2_temporal.o
	@mkdir -p $(@D)
	$(FC) $(FFLAGS_CHECK) -I $(BUILD_DIR)/interfaces $^ -o $@

m9d3-spatial-check: $(BUILD_DIR)/m9d3_spatial_regression

$(BUILD_DIR)/m9d3_spatial_regression: $(BUILD_DIR)/obj/check/src/fortran77/models/eos_constant_z.o $(BUILD_DIR)/obj/check/src/fortran77/models/friction_swamee_jain.o $(BUILD_DIR)/obj/check/src/fortran77/transient_residual.o $(BUILD_DIR)/obj/check/src/fortran77/transient_jacobian.o $(BUILD_DIR)/obj/check/src/fortran77/linear_solve.o $(BUILD_DIR)/obj/check/src/fortran77/input_validation.o $(BUILD_DIR)/obj/check/src/fortran77/newton_solver.o $(BUILD_DIR)/obj/check/src/fortran77/transient_step.o $(BUILD_DIR)/obj/check/src/fortran77/transient_integrate.o $(BUILD_DIR)/obj/check/tests/fortran77/m9d3_spatial_regression.o
	@mkdir -p $(@D)
	$(FC) $(FFLAGS_CHECK) -I $(BUILD_DIR)/interfaces $^ -o $@

m9d4-inventory-check: $(BUILD_DIR)/m9d4_inventory_regression

$(BUILD_DIR)/m9d4_inventory_regression: $(BUILD_DIR)/obj/check/src/fortran77/models/eos_constant_z.o $(BUILD_DIR)/obj/check/src/fortran77/models/friction_swamee_jain.o $(BUILD_DIR)/obj/check/src/fortran77/transient_residual.o $(BUILD_DIR)/obj/check/src/fortran77/transient_jacobian.o $(BUILD_DIR)/obj/check/src/fortran77/linear_solve.o $(BUILD_DIR)/obj/check/src/fortran77/input_validation.o $(BUILD_DIR)/obj/check/src/fortran77/newton_solver.o $(BUILD_DIR)/obj/check/src/fortran77/transient_step.o $(BUILD_DIR)/obj/check/src/fortran77/transient_integrate.o $(BUILD_DIR)/obj/check/tests/fortran77/m9d4_inventory_regression.o
	@mkdir -p $(@D)
	$(FC) $(FFLAGS_CHECK) -I $(BUILD_DIR)/interfaces $^ -o $@

m9d5-halfcell-check: $(BUILD_DIR)/m9d5_halfcell_regression

$(BUILD_DIR)/m9d5_halfcell_regression: $(BUILD_DIR)/obj/check/src/fortran77/models/eos_constant_z.o $(BUILD_DIR)/obj/check/src/fortran77/models/friction_swamee_jain.o $(BUILD_DIR)/obj/check/src/fortran77/transient_residual.o $(BUILD_DIR)/obj/check/src/fortran77/transient_jacobian.o $(BUILD_DIR)/obj/check/tests/fortran77/m9d5_halfcell_regression.o
	@mkdir -p $(@D)
	$(FC) $(FFLAGS_CHECK) -I $(BUILD_DIR)/interfaces $^ -o $@

m9e2-interface-check: $(BUILD_DIR)/m9e2_interface_check

$(BUILD_DIR)/m9e2_interface_check: $(BUILD_DIR)/obj/check/src/interfaces/pipe_solver_api.o $(BUILD_DIR)/obj/check/src/fortran77/models/eos_constant_z.o $(BUILD_DIR)/obj/check/src/fortran77/models/friction_swamee_jain.o $(BUILD_DIR)/obj/check/src/fortran77/transient_residual.o $(BUILD_DIR)/obj/check/src/fortran77/transient_jacobian.o $(BUILD_DIR)/obj/check/src/fortran77/linear_solve.o $(BUILD_DIR)/obj/check/src/fortran77/input_validation.o $(BUILD_DIR)/obj/check/src/fortran77/newton_solver.o $(BUILD_DIR)/obj/check/src/fortran77/transient_step.o $(BUILD_DIR)/obj/check/src/fortran77/transient_integrate.o $(BUILD_DIR)/obj/check/src/fortran77/transient_adaptive.o $(BUILD_DIR)/obj/check/src/fortran77/transient_adaptive_schedule.o $(BUILD_DIR)/obj/check/tests/fortran77/m9e2_interface_check.o
	@mkdir -p $(@D)
	$(FC) $(FFLAGS_CHECK) -I $(BUILD_DIR)/interfaces $^ -o $@

m9e2b-integration-check: $(BUILD_DIR)/m9e2b_integration_interface

$(BUILD_DIR)/m9e2b_integration_interface: $(BUILD_DIR)/obj/check/src/interfaces/pipe_solver_api.o $(BUILD_DIR)/obj/check/src/fortran77/models/eos_constant_z.o $(BUILD_DIR)/obj/check/src/fortran77/models/friction_swamee_jain.o $(BUILD_DIR)/obj/check/src/fortran77/transient_residual.o $(BUILD_DIR)/obj/check/src/fortran77/transient_jacobian.o $(BUILD_DIR)/obj/check/src/fortran77/linear_solve.o $(BUILD_DIR)/obj/check/src/fortran77/input_validation.o $(BUILD_DIR)/obj/check/src/fortran77/newton_solver.o $(BUILD_DIR)/obj/check/src/fortran77/transient_step.o $(BUILD_DIR)/obj/check/src/fortran77/transient_integrate.o $(BUILD_DIR)/obj/check/src/fortran77/transient_adaptive.o $(BUILD_DIR)/obj/check/src/fortran77/transient_adaptive_schedule.o $(BUILD_DIR)/obj/check/tests/fortran77/m9e2b_integration_interface.o
	@mkdir -p $(@D)
	$(FC) $(FFLAGS_CHECK) -I $(BUILD_DIR)/interfaces $^ -o $@

m9e2c-adaptive-check: $(BUILD_DIR)/m9e2c_adaptive_interface

$(BUILD_DIR)/m9e2c_adaptive_interface: $(BUILD_DIR)/obj/check/src/interfaces/pipe_solver_api.o $(BUILD_DIR)/obj/check/src/fortran77/models/eos_constant_z.o $(BUILD_DIR)/obj/check/src/fortran77/models/friction_swamee_jain.o $(BUILD_DIR)/obj/check/src/fortran77/transient_residual.o $(BUILD_DIR)/obj/check/src/fortran77/transient_jacobian.o $(BUILD_DIR)/obj/check/src/fortran77/linear_solve.o $(BUILD_DIR)/obj/check/src/fortran77/input_validation.o $(BUILD_DIR)/obj/check/src/fortran77/newton_solver.o $(BUILD_DIR)/obj/check/src/fortran77/transient_step.o $(BUILD_DIR)/obj/check/src/fortran77/transient_integrate.o $(BUILD_DIR)/obj/check/src/fortran77/transient_adaptive.o $(BUILD_DIR)/obj/check/src/fortran77/transient_adaptive_schedule.o $(BUILD_DIR)/obj/check/tests/fortran77/m9e2c_adaptive_interface.o
	@mkdir -p $(@D)
	$(FC) $(FFLAGS_CHECK) -I $(BUILD_DIR)/interfaces $^ -o $@

m9e2d-schedule-check: $(BUILD_DIR)/m9e2d_schedule_interface

$(BUILD_DIR)/m9e2d_schedule_interface: $(BUILD_DIR)/obj/check/src/interfaces/pipe_solver_api.o $(BUILD_DIR)/obj/check/src/fortran77/models/eos_constant_z.o $(BUILD_DIR)/obj/check/src/fortran77/models/friction_swamee_jain.o $(BUILD_DIR)/obj/check/src/fortran77/transient_residual.o $(BUILD_DIR)/obj/check/src/fortran77/transient_jacobian.o $(BUILD_DIR)/obj/check/src/fortran77/linear_solve.o $(BUILD_DIR)/obj/check/src/fortran77/input_validation.o $(BUILD_DIR)/obj/check/src/fortran77/newton_solver.o $(BUILD_DIR)/obj/check/src/fortran77/transient_step.o $(BUILD_DIR)/obj/check/src/fortran77/transient_integrate.o $(BUILD_DIR)/obj/check/src/fortran77/transient_adaptive.o $(BUILD_DIR)/obj/check/src/fortran77/transient_adaptive_schedule.o $(BUILD_DIR)/obj/check/tests/fortran77/m9e2d_schedule_interface.o
	@mkdir -p $(@D)
	$(FC) $(FFLAGS_CHECK) -I $(BUILD_DIR)/interfaces $^ -o $@

check: m9e3-failure-check steady-residual friction-check jacobian-check newton-check timestep-check integration-check adaptive-check adaptive-schedule-check adaptive-schedule-pressure-check input-validation-check m9b3-negative-check m9c2a-failure-check m9c2b-failure-check m9c3-partial-failure-check m9d2-temporal-check m9d3-spatial-check m9d4-inventory-check m9d5-halfcell-check m9e2-interface-check m9e2b-integration-check m9e2c-adaptive-check m9e2d-schedule-check
	./$(BUILD_DIR)/m9e3_failure_interface
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
	./$(BUILD_DIR)/m9e2b_integration_interface
	./$(BUILD_DIR)/m9e2c_adaptive_interface
	./$(BUILD_DIR)/m9e2d_schedule_interface

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

m9c3-partial-failure-check: $(BUILD_DIR)/m9c3_partial_failure_check

$(BUILD_DIR)/m9c3_partial_failure_check: $(BUILD_DIR)/obj/check/src/fortran77/models/eos_constant_z.o $(BUILD_DIR)/obj/check/src/fortran77/models/friction_swamee_jain.o $(BUILD_DIR)/obj/check/src/fortran77/transient_residual.o $(BUILD_DIR)/obj/check/src/fortran77/input_validation.o $(BUILD_DIR)/obj/check/src/fortran77/transient_integrate.o $(BUILD_DIR)/obj/check/src/fortran77/transient_adaptive.o $(BUILD_DIR)/obj/check/src/fortran77/transient_adaptive_schedule.o $(BUILD_DIR)/obj/check/tests/fortran77/m9c3_partial_failure_check.o
	@mkdir -p $(@D)
	$(FC) $(FFLAGS_CHECK) -I $(BUILD_DIR)/interfaces $^ -o $@


# M9E.4: debug objects are reusable across regression executables.
$(BUILD_DIR)/obj/check/%.o: %.f
	@mkdir -p $(@D)
	$(FC) $(FFLAGS_CHECK) -c $< -o $@

$(BUILD_DIR)/obj/check/%.o: %.f90
	@mkdir -p $(@D)
	$(FC) $(FFLAGS_CHECK) -I $(BUILD_DIR)/interfaces -c $< -o $@

$(BUILD_DIR)/obj/check/src/interfaces/pipe_solver_api.o: src/interfaces/pipe_solver_api.f90
	@mkdir -p $(@D) $(BUILD_DIR)/interfaces
	$(FC) $(FFLAGS_CHECK) -J $(BUILD_DIR)/interfaces -c $< -o $@

$(BUILD_DIR)/obj/check/tests/fortran77/m9e2_interface_check.o: $(BUILD_DIR)/obj/check/src/interfaces/pipe_solver_api.o
$(BUILD_DIR)/obj/check/tests/fortran77/m9e2b_integration_interface.o: $(BUILD_DIR)/obj/check/src/interfaces/pipe_solver_api.o
$(BUILD_DIR)/obj/check/tests/fortran77/m9e2c_adaptive_interface.o: $(BUILD_DIR)/obj/check/src/interfaces/pipe_solver_api.o
$(BUILD_DIR)/obj/check/tests/fortran77/m9e2d_schedule_interface.o: $(BUILD_DIR)/obj/check/src/interfaces/pipe_solver_api.o

$(BUILD_DIR)/obj/check/tests/fortran77/m9e3_failure_interface.o: $(BUILD_DIR)/obj/check/src/interfaces/pipe_solver_api.o

m9e3-failure-check: $(BUILD_DIR)/m9e3_failure_interface

$(BUILD_DIR)/m9e3_failure_interface: $(BUILD_DIR)/obj/check/src/interfaces/pipe_solver_api.o $(BUILD_DIR)/obj/check/src/fortran77/models/eos_constant_z.o $(BUILD_DIR)/obj/check/src/fortran77/models/friction_swamee_jain.o $(BUILD_DIR)/obj/check/src/fortran77/transient_residual.o $(BUILD_DIR)/obj/check/src/fortran77/transient_jacobian.o $(BUILD_DIR)/obj/check/src/fortran77/linear_solve.o $(BUILD_DIR)/obj/check/src/fortran77/input_validation.o $(BUILD_DIR)/obj/check/src/fortran77/newton_solver.o $(BUILD_DIR)/obj/check/src/fortran77/transient_step.o $(BUILD_DIR)/obj/check/src/fortran77/transient_integrate.o $(BUILD_DIR)/obj/check/src/fortran77/transient_adaptive.o $(BUILD_DIR)/obj/check/src/fortran77/transient_adaptive_schedule.o $(BUILD_DIR)/obj/check/tests/fortran77/m9e3_failure_interface.o
	$(FC) $(FFLAGS_CHECK) -I $(BUILD_DIR)/interfaces $^ -o $@

include m10.mk
