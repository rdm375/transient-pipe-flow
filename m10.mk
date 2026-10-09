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

# M10.2: component profiling (no changes to production Fortran sources).
M10_PROFILE_CORE = $(filter-out src/fortran77/newton_solver.f,$(M10_CORE))

.PHONY: m10-profile-build

m10-profile-build: build/m10_cpu_profile

build/m10_cpu_profile: $(M10_CORE) benchmarks/m10_cpu_baseline.f90 \
    benchmarks/m10_profile.f90 benchmarks/generate_m10_profile.py
	@mkdir -p build
	python3 benchmarks/generate_m10_profile.py
	$(FC) $(FFLAGS_OPT) -Wall -Wextra \
	    benchmarks/m10_profile.f90 $(M10_PROFILE_CORE) \
	    build/m10_newton_profile.f build/m10_cpu_profile.f90 -o $@
