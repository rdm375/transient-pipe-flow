#pragma once

namespace pipe_sim {

struct FrictionValueDerivative {
    double value;
    double derivative_re;
};

double density_cz(double pressure, double temperature,
                  double z, double rs);

double density_derivative_cz(double temperature,
                             double z, double rs);

double reynolds_mass(double mdot, double diameter,
                     double area, double viscosity);

double reynolds_derivative_mdot(double mdot, double diameter,
                               double area, double viscosity);

double friction_sj(double reynolds, double roughness,
                   double diameter);

double friction_derivative_re(double reynolds, double roughness,
                              double diameter);

FrictionValueDerivative friction_sj_value_derivative(
    double reynolds, double roughness, double diameter);

double friction_source(double fd, double diameter,
                       double area, double density_face,
                       double mdot);

double friction_source_derivative_mdot(
    double fd, double dfdm, double diameter,
    double area, double density_face, double mdot);

} // namespace pipe_sim
