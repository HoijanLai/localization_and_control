#include "MPC.h"
#include "utils.h"
#include <iostream>
#include <fstream>
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#include "json.hpp"

// for convenience
using json = nlohmann::json;
using CppAD::AD;




// Set the timestep length and duration
size_t N = 11;
double dt = 0.1;

// set the vehicle parameters
const double Lf = 2.67;

// speed reference against stop
const double v_ref = 40;

// index for convenience
size_t x_start     = 0;
size_t y_start     = x_start + N;
size_t psi_start   = y_start + N;
size_t v_start     = psi_start + N;
size_t cte_start   = v_start + N;
size_t epsi_start  = cte_start + N;
size_t delta_start = epsi_start + N;
size_t a_start     = delta_start + N - 1;

// Weights tuned for the cost
double w_cte   = 1.0;
double w_epsi  = 10.0;
double w_v     = 1.0;
double w_delta = 200.0;
double w_a     = 1.0;
double w_dd    = 20.0;
double w_da    = 40.0;


/* ===========================================

 The fancy definition of this FG_eval system

=========================================== */

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;


  void operator()(ADvector& fg, const ADvector& vars) {

    /* ----------------
      Cost Definition
    ---------------- */

    // Minimize the error
    fg[0] = 0;
    for (unsigned int t = 0; t < N; ++t) {
      fg[0] += w_cte * CppAD::pow(vars[cte_start + t], 2);
      fg[0] += w_epsi * CppAD::pow(vars[epsi_start + t], 2);
      fg[0] += w_v * CppAD::pow(vars[v_start + t] - v_ref, 2);
    }

    // Minimize the use of actuators.
    for (unsigned int t = 0; t < N - 1; t++) {
      fg[0] += w_delta * CppAD::pow(vars[delta_start + t], 2);
      fg[0] += w_a * CppAD::pow(vars[a_start + t], 2);
    }

    // Minimize the value gap between sequential actuations.
    for (unsigned int t = 0; t < N - 2; t++) {
      fg[0] += w_dd * CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
      fg[0] += w_da * CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
    }

    /* ------------
      Constraints
    ------------ */

    // The initial constraints
    fg[1 + x_start]    = vars[x_start];
    fg[1 + y_start]    = vars[y_start];
    fg[1 + psi_start]  = vars[psi_start];
    fg[1 + v_start]    = vars[v_start];
    fg[1 + cte_start]  = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    // The rest of the constraints
    for (unsigned int t = 1; t < N; t++) {
      // The state at time t+1 .
      AD<double> x1    = vars[x_start + t];
      AD<double> y1    = vars[y_start + t];
      AD<double> psi1  = vars[psi_start + t];
      AD<double> v1    = vars[v_start + t];
      AD<double> cte1  = vars[cte_start + t];
      AD<double> epsi1 = vars[epsi_start + t];

      // The state at time t.
      AD<double> x0    = vars[x_start + t - 1];
      AD<double> y0    = vars[y_start + t - 1];
      AD<double> psi0  = vars[psi_start + t - 1];
      AD<double> v0    = vars[v_start + t - 1];
      AD<double> cte0  = vars[cte_start + t - 1];
      AD<double> epsi0 = vars[epsi_start + t - 1];

      // Only consider the actuation at time t.
      AD<double> delta0 = vars[delta_start + t - 1];
      AD<double> a0 = vars[a_start + t - 1];

      // The desired state vars for y and psi
      AD<double> f0 = polyeval(coeffs, x0);
      AD<double> psides0 = CppAD::atan(polyderivative(coeffs, x0));

      // the constraints are set based on vehicle motion
      fg[1 + x_start + t]   = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_start + t]   = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_start + t] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
      fg[1 + v_start + t]   = v1 - (v0 + a0 * dt);

      fg[1 + cte_start + t] =
      cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));

      fg[1 + epsi_start + t] =
      epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
    }
  }
};




/* =====================================

   MPC class definition implementation

===================================== */

MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {

  /* ----------------------
    Just some preparation
  ---------------------- */

  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // Set the number of model variables (includes both states and inputs).
  // & the number of constraints
  size_t n_vars = N * 6 + 2 * (N - 1);
  size_t n_constraints = N * 6;


  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (unsigned int i = 0; i < n_vars; i++) {
    vars[i] = 0.0;
  }
  vars[x_start]    = state[0];
  vars[y_start]    = state[1];
  vars[psi_start]  = state[2];
  vars[v_start]    = state[3];
  vars[cte_start]  = state[4];
  vars[epsi_start] = state[5];

  /* ---------------------------------------------
    Set the lower and upper bounds for variables
  --------------------------------------------- */

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (unsigned int i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  for (unsigned int i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }

  // Acceleration/decceleration upper and lower limits.
  for (unsigned int i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }



  /* ------------------------------------------------
     Set the lower and upper bounds for constraints
     0 besides initial state.
  ------------------------------------------------ */

  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (unsigned int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  constraints_lowerbound[x_start]    = state[0];
  constraints_lowerbound[y_start]    = state[1];
  constraints_lowerbound[psi_start]  = state[2];
  constraints_lowerbound[v_start]    = state[3];
  constraints_lowerbound[cte_start]  = state[4];
  constraints_lowerbound[epsi_start] = state[5];

  constraints_upperbound[x_start]    = state[0];
  constraints_upperbound[y_start]    = state[1];
  constraints_upperbound[psi_start]  = state[2];
  constraints_upperbound[v_start]    = state[3];
  constraints_upperbound[cte_start]  = state[4];
  constraints_upperbound[epsi_start] = state[5];


  /* -------------------------------------------------
    Solve the optimization and return the actuations
  ------------------------------------------------- */

  // Object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  // Options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  options += "Numeric max_cpu_time          0.5\n";

  // Place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // Solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;
  // Return the first two actuations
  std::vector<double> sol = {solution.x[delta_start], solution.x[a_start]};
  for(unsigned int t = 0; t < N; t++) {
    sol.push_back(solution.x[x_start + t]);
    sol.push_back(solution.x[y_start + t]);
  }
  return sol;

}
// End MPC::Solve
