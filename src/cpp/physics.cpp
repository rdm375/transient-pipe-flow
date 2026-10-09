#include "pipe_sim/physics.hpp"

#include <cmath>
#include <stdexcept>

namespace pipe_sim {

double density_cz(double pressure, double temperature,
                  double z, double rs)
{
    return pressure / (z * rs * temperature);
}

double density_derivative_cz(double temperature,
                             double z, double rs)
{
    return 1.0 / (z * rs * temperature);
}

double reynolds_mass(double mdot, double diameter,
                     double area, double viscosity)
{
    if (!(diameter > 0.0 && area > 0.0 && viscosity > 0.0))
        throw std::invalid_argument("Invalid Reynolds parameters");

    return std::abs(mdot) * diameter / (area * viscosity);
}

double reynolds_derivative_mdot(double mdot, double diameter,
                               double area, double viscosity)
{
    if (!(diameter > 0.0 && area > 0.0 && viscosity > 0.0))
        throw std::invalid_argument("Invalid Reynolds parameters");

    return std::copysign(diameter / (area * viscosity), mdot);
}

double friction_sj(double reynolds, double roughness,
                   double diameter)
{
    if (!(reynolds > 0.0 && roughness >= 0.0 && diameter > 0.0))
        throw std::invalid_argument("Invalid friction parameters");

    const double arg =
        roughness / (3.7 * diameter)
        + 5.74 / std::pow(reynolds, 0.9);

    const double larg = std::log10(arg);

    return 0.25 / (larg * larg);
}

double friction_derivative_re(double reynolds, double roughness,
                              double diameter)
{
    if (!(reynolds > 0.0 && roughness >= 0.0 && diameter > 0.0))
        throw std::invalid_argument("Invalid friction parameters");

    const double arg =
        roughness / (3.7 * diameter)
        + 5.74 / std::pow(reynolds, 0.9);

    const double larg = std::log10(arg);
    const double ln10 = std::log(10.0);
    const double dxdre =
        -0.9 * 5.74 / std::pow(reynolds, 1.9);

    return -0.5 * dxdre /
           (arg * ln10 * larg * larg * larg);
}

FrictionValueDerivative friction_sj_value_derivative(
    double reynolds, double roughness, double diameter)
{
    if (!(reynolds > 0.0 && roughness >= 0.0 && diameter > 0.0))
        throw std::invalid_argument("Invalid friction parameters");

    const double re09 = std::pow(reynolds, 0.9);
    const double arg =
        roughness / (3.7 * diameter) + 5.74 / re09;

    const double larg = std::log10(arg);
    const double fd = 0.25 / (larg * larg);

    const double ln10 = std::log(10.0);
    const double dxdre = -0.9 * 5.74 / (reynolds * re09);

    const double dfdre =
        -0.5 * dxdre / (arg * ln10 * larg * larg * larg);

    return {fd, dfdre};
}

double friction_source(double fd, double diameter,
                       double area, double density_face,
                       double mdot)
{
    if (!(diameter > 0.0 && area > 0.0 && density_face > 0.0))
        throw std::invalid_argument("Invalid friction source parameters");

    return fd * mdot * std::abs(mdot) /
           (2.0 * diameter * area * area * density_face);
}

double friction_source_derivative_mdot(
    double fd, double dfdm, double diameter,
    double area, double density_face, double mdot)
{
    if (!(diameter > 0.0 && area > 0.0 && density_face > 0.0))
        throw std::invalid_argument("Invalid friction source parameters");

    const double c =
        1.0 / (2.0 * diameter * area * area * density_face);

    return c * (dfdm * mdot * std::abs(mdot)
                + 2.0 * fd * std::abs(mdot));
}

} // namespace pipe_sim
