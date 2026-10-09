#include "pipe_sim/integrate.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <vector>

using namespace pipe_sim;

int main(int argc, char** argv)
{
    const int steps = argc > 1 ? std::atoi(argv[1]) : 80;
    const int batch = argc > 2 ? std::atoi(argv[2]) : 2048;

    if (steps < 0 || batch < 1) {
        std::fprintf(stderr, "Invalid steps or batch size\n");
        return 1;
    }

    constexpr int n = 100;
    constexpr double inlet = 100.0;

    const double pi = std::acos(-1.0);
    const double diameter = 1.0;
    const double area = pi * diameter * diameter / 4.0;

    const TransientParameters par{
        100000.0 / n, 15.0, 0.65,
        diameter, area, 288.15, 0.9,
        500.0, 1.1e-5, 4.5e-5
    };

    const DemandRamp demand{
        8.0e6, 100.0, 105.0, 600.0
    };

    const NewtonOptions opt{
        1.0e-9, 1.0e-12, 30, false
    };

    std::vector<double> initial_pressure(n + 1);
    std::vector<double> initial_flow(n, inlet);

    const double reynolds =
        inlet * diameter / (area * par.viscosity);

    const double friction_arg =
        par.roughness / (3.7 * diameter)
        + 5.74 / std::pow(reynolds, 0.9);

    const double friction =
        0.25 / std::pow(std::log10(friction_arg), 2.0);

    const double coefficient =
        friction * par.z * par.gas_constant
        * par.temperature * inlet * inlet
        / (diameter * area * area);

    for (int i = 0; i <= n; ++i) {
        initial_pressure[i] = std::sqrt(
            demand.inlet_pressure * demand.inlet_pressure
            - coefficient * i * par.dx
        );
    }

    // Establish an independent numerical reference.
    auto reference_pressure = initial_pressure;
    auto reference_flow = initial_flow;
    double reference_inlet = inlet;

    std::vector<double> rt(steps + 1), ri(steps + 1);
    std::vector<double> ro(steps + 1), rp(steps + 1);
    std::vector<double> rl(steps + 1);
    std::vector<int> rn(steps + 1);

    IntegrationHistory reference_history{
        rt, ri, ro, rp, rl, rn
    };

    IntegrationWorkspace reference_workspace(n);

    const auto reference = integrate_transient(
        n, reference_pressure, reference_flow,
        reference_inlet, demand, steps,
        par, opt, reference_history, reference_workspace
    );

    if (reference.info != 0 ||
        reference.completed_steps != steps) {
        std::fprintf(stderr, "CPU reference failed\n");
        return 1;
    }

    for (int workers : {1, 2, 4, 8}) {
        std::atomic<int> failures{0};
        std::atomic<int> next_pipe{0};

        // Construct private workspaces before timing.
        struct Worker {
            IntegrationWorkspace workspace;
            std::vector<double> pressure;
            std::vector<double> flow;
            std::vector<double> time, hin, hout;
            std::vector<double> hp, hl;
            std::vector<int> iterations;

            Worker(int n, int steps)
                : workspace(n),
                  pressure(n + 1),
                  flow(n),
                  time(steps + 1),
                  hin(steps + 1),
                  hout(steps + 1),
                  hp(steps + 1),
                  hl(steps + 1),
                  iterations(steps + 1)
            {}
        };

        std::vector<Worker> contexts;
        contexts.reserve(workers);

        for (int w = 0; w < workers; ++w)
            contexts.emplace_back(n, steps);

        auto execute_batch = [&](bool validate) {
            next_pipe.store(0, std::memory_order_relaxed);

        std::vector<std::thread> threads;
        threads.reserve(workers);

        for (int w = 0; w < workers; ++w) {
            threads.emplace_back([&, w]() {
                auto& ctx = contexts[w];

                const IntegrationHistory history{
                    ctx.time, ctx.hin, ctx.hout,
                    ctx.hp, ctx.hl, ctx.iterations
                };

                while (true) {
                    const int pipe = next_pipe.fetch_add(
                        1, std::memory_order_relaxed
                    );

                    if (pipe >= batch)
                        break;

                    std::copy(
                        initial_pressure.begin(),
                        initial_pressure.end(),
                        ctx.pressure.begin()
                    );

                    std::copy(
                        initial_flow.begin(),
                        initial_flow.end(),
                        ctx.flow.begin()
                    );

                    double current_inlet = inlet;

                    const auto result = integrate_transient(
                        n, ctx.pressure, ctx.flow,
                        current_inlet, demand, steps,
                        par, opt, history, ctx.workspace
                    );

                    if (validate) {
                    if (result.info != reference.info ||
                        result.completed_steps !=
                            reference.completed_steps) {
                        failures.fetch_add(
                            1, std::memory_order_relaxed
                        );
                        continue;
                    }

                    // Exact CPU-to-CPU validation.
                    bool mismatch =
                        current_inlet != reference_inlet ||
                        ctx.pressure != reference_pressure ||
                        ctx.flow != reference_flow ||
                        result.max_mass_defect !=
                            reference.max_mass_defect ||
                        ctx.iterations != rn ||
                        ctx.hin != ri ||
                        ctx.hout != ro ||
                        ctx.hp != rp ||
                        ctx.hl != rl ||
                        ctx.time != rt;

                    if (mismatch)
                        failures.fetch_add(
                            1, std::memory_order_relaxed
                        );
                    }

                }
            });
        }

        for (auto& thread : threads)
            thread.join();
        };

        const auto start = std::chrono::steady_clock::now();

        execute_batch(false);

        const auto stop = std::chrono::steady_clock::now();

        // Validation is deliberately outside the timed interval.
        execute_batch(true);

        const double ms =
            std::chrono::duration<double, std::milli>(
                stop - start
            ).count();

        std::printf(
            "CPU_BATCH workers=%d batch=%d steps=%d "
            "elapsed_ms=%.4f simulations_per_second=%.3f "
            "failures=%d\n",
            workers, batch, steps, ms,
            1000.0 * batch / ms,
            failures.load()
        );

        if (failures.load() != 0)
            return 1;
    }

    std::puts("PASS M14.6c multicore CPU batch");
    return 0;
}
