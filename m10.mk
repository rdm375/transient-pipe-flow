# M10.1: optimized end-to-end CPU baseline.
M10_CORE = \
    src/fortran77/models/eos_constant_z.f \
    src/fortran77/models/friction_swamee_jain.f \
    src/fortran77/transient_residual.f \
    src/fortran77/transient_jacobian.f \
    src/fortran77/linear_solve.f \
    src/fortran77/input_validation.f \
    src/fortran77/newton_solver.f \
    src/fortran77/transient_step.f \
    src/fortran77/transient_integrate.f

.PHONY: m10-baseline-build m10-baseline

m10-baseline-build: build/m10_cpu_baseline

build/m10_cpu_baseline: $(M10_CORE) benchmarks/m10_cpu_baseline.f90
	@mkdir -p build
	$(FC) $(FFLAGS_OPT) -Wall -Wextra $(M10_CORE) \
	    benchmarks/m10_cpu_baseline.f90 -o $@

m10-baseline: m10-baseline-build
	python3 benchmarks/run_m10_baseline.py
