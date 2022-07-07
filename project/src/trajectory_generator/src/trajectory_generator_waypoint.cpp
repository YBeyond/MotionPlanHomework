#include "trajectory_generator_waypoint.h"
#include "osqp++.h"
#include <fstream>
#include <iostream>
#include <ros/console.h>
#include <ros/ros.h>
#include <stdio.h>
#include <string>

using namespace std;
using namespace Eigen;

#define inf 1 >> 30

TrajectoryGeneratorWaypoint::TrajectoryGeneratorWaypoint() {}

TrajectoryGeneratorWaypoint::~TrajectoryGeneratorWaypoint() {}

/**
  * @brief get partial factorial generated in derivative computing
  *
  * @param[in] N polynomial term power
  * @param[in] D derivative order
  *
  * @return partial factorial as N!/D!
  */
double TrajectoryGeneratorWaypoint::GetFactorial(const int N, const int D) {
    double result{1.0};

    for (int i = 0; i < D; ++i) {
        result *= (N - i);
    }

    return result;
}

/**
 * @brief generate minimum snap trajectory through numeric method with OSQP C++
 *
 * @param[in] tOrder the L2-norm of [tOrder]th derivative of the target trajectory will be used as objective function
 * @param[in] cOrder continuity constraints, the target trajectory should be [cOrder]th continuous at intermediate waypoints
 * @param[in] Pos equality constraints, the target trajectory should pass all the waypoints, (K + 1)-by-1
 * @param[in] Vel equality constraints, boundary(start & goal) velocity specifications, 2-by-1
 * @param[in] Acc equality constraints, boundary(start & goal) acceleration specifications, 2-by-1
 * @param[in] Time pre-computed time allocations, K-by-1
 *
 * @return polynomial coeffs of generated trajectory, K-by-N
 * @note the pre-assumption is no allocated segment time in Time is 0
 */
Eigen::MatrixXd TrajectoryGeneratorWaypoint::PolyQPGenerationNumeric(
    const int tOrder,       
    const int cOrder,    
    const Eigen::VectorXd &Pos,
    const Eigen::VectorXd &Vel,
    const Eigen::VectorXd &Acc,
    const Eigen::VectorXd &Time
) {
    // num. of polynomial coeffs:
    const auto N = (cOrder + 1) << 1;  // [cOrder]阶连续，取值为3时，为7阶多项式，8个系数
    // num. of trajectory segments:
    const auto K = Time.size();
    // dim of flattened output:
    const auto D = K * N;
    // num. of inequality constraints:
    const auto C = (
        // 1. boundary value equality constraints:
        N +
        // 2. intermediate waypoint passing equality constraints:
        (K - 1) +
        // 3. intermediate waypoint continuity constaints:
        (K - 1) * (cOrder + 1));
    // 点的维度（x,y,z）

    //
    // init output:
    //
    //     the trajectory segment k is defined by result(k, 0) + t*(result(k, 1) + ...)
    // 
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(K, N);

    //
    // problem definition:
    // 
    //     min: 0.5 * x'Px + q'x
    //     s.t. l <= Ax <= u
    // 
    // osqp-cpp interface
    //
    //     P: objective_matrix
    //     q: objective_vector
    //     A: constraint_matrix
    //     l: lower_bounds
    //     u: upper_bounds
    // 

    // 
    // define objective matrix P:
    //
    // 1. init:
    Eigen::SparseMatrix<double> P(D, D);
    std::vector<Eigen::Triplet<double>> PTriplets;
    
    // 2. cache results for objective matrix construction:
    
    // 2.a time power
    Eigen::VectorXd PTimePower = Eigen::VectorXd::Ones(K);
    for (int c = 1; c < tOrder; ++c) {
        PTimePower = PTimePower.cwiseProduct(Time);
    }
    
    // 2.b factorial
    std::map<int, double> PFactorial;
    for (int n = tOrder; n < N; ++n) {
        PFactorial.insert(
            {n, GetFactorial(n, tOrder)}
        );
    }

    // 3. populate PTriplets:
    for (int k = 0; k < K; ++k) {
        const auto currentSegmentIdxOffset = k * N;

        for (int m = tOrder; m < N; ++m) {
            for (int n = tOrder; n < N; ++n) {
                    PTriplets.emplace_back(
                    currentSegmentIdxOffset + m,
                    currentSegmentIdxOffset + n,
                    Time(k)*PFactorial[m]*PFactorial[n]/(
                        (m + n - (tOrder << 1) + 1) *
                        PTimePower(k) * PTimePower(k)        // t这一项是有问题的？
                    )
                );
            }
        }
    }

    // 4. populate P:
    P.setFromTriplets(std::begin(PTriplets),std::end(PTriplets));

    //
    // define constraint matrix A:
    //
    // 1. init:
    Eigen::SparseMatrix<double> A(C, D);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(C);
    std::vector<Eigen::Triplet<double>> ATriplets;

    // 2. cache results for constraint matrix construction:
    std::map<int, Eigen::VectorXd> ATimePower;
    std::map<int, double> AFactorial;

    // 2.a init:
    ATimePower.insert(
        {0, Eigen::VectorXd::Ones(K)}
    );
    
    auto GetAFactorialKey = [&](const int n, const int d) {
        return d * N + n;
    };
    for (int n = 0; n < N; ++n) {
        AFactorial.insert(
            {GetAFactorialKey(n, 0), 1.0}
        );
    }

    for (int c = 1; c <= cOrder; ++c) {
        // 2.b time power
        ATimePower.insert(
            {c, ATimePower[c - 1].cwiseProduct(Time)}
        );

        // 2.c factorial
        for (int n = c; n < N; ++n) {
            AFactorial.insert(
                {GetAFactorialKey(n, c), (n - c + 1)*AFactorial[GetAFactorialKey(n, c - 1)]}
            );
        }
    };

    // 3. populate ATriplets:
    {
        int currentConstraintIdx{0};

        //
        // 3.1 boundary value equality constraints:
        //

        // init:
        Eigen::MatrixXd boundaryValues = Eigen::MatrixXd::Zero(2, cOrder + 1);
        boundaryValues(0, 0) = Pos(0); // 两行三列？点不应该是一个3D点吗？
        boundaryValues(0, 1) = Vel(0);
        boundaryValues(0, 2) = Acc(0);
        boundaryValues(1, 0) = Pos(K);
        boundaryValues(1, 1) = Vel(1);
        boundaryValues(1, 2) = Acc(1);

        // populate constraints:
        for (int c = 0; c <= cOrder; ++c) {
            // start waypoint:
            ATriplets.emplace_back(
                currentConstraintIdx, c, AFactorial[GetAFactorialKey(c, c)]/ATimePower[c](0)
            );
            b(currentConstraintIdx) = boundaryValues(0, c);
            
            ++currentConstraintIdx;

            // end waypoint:
            for (int n = c; n < N; ++n) {
                ATriplets.emplace_back(
                    currentConstraintIdx, (K - 1)*N + n, AFactorial[GetAFactorialKey(n, c)]/ATimePower[c](K - 1)
                );
            }
            b(currentConstraintIdx) = boundaryValues(1, c);

            ++currentConstraintIdx;
        }

        //
        // 3.2 intermediate waypoint passing equality constraints:
        //
        for (int k = 1; k < K; ++k) {
            ATriplets.emplace_back(
                currentConstraintIdx, k*N, 1.0
            );
            b(currentConstraintIdx) = Pos(k);

            ++currentConstraintIdx;
        }

        //
        // 3.3 intermediate waypoint continuity constaints:
        //
        for (int c = 0; c <= cOrder; ++c) {
            for (int k = 1; k < K; ++k) {
                // the end of previous trajectory segment:
                for (int n = c; n < N; ++n) {
                    ATriplets.emplace_back(
                        currentConstraintIdx, 
                        (k - 1) * N + n,
                        AFactorial[GetAFactorialKey(n, c)]/ATimePower[c](k - 1)
                    );
                }

                // should equal to the start of current trajectory segment:
                ATriplets.emplace_back(
                    currentConstraintIdx,
                    k * N + c,
                    -AFactorial[GetAFactorialKey(c, c)]/ATimePower[c](k) 
                );

                // set next constraint:
                ++currentConstraintIdx;
            }
        }
    }

    // 4. populate A:
    A.setFromTriplets(std::begin(ATriplets),std::end(ATriplets));

    //
    // solve with OSQP c++
    //
    osqp::OsqpInstance instance;
    instance.objective_matrix = P;
    instance.objective_vector = Eigen::VectorXd::Zero(D);

    instance.constraint_matrix = A;
    instance.lower_bounds = instance.upper_bounds = b;

    osqp::OsqpSolver solver;
    osqp::OsqpSettings settings;

    // init solver:
    if (solver.Init(instance, settings).ok()) {
        // solve.
        const auto exitCode = solver.Solve();

        if (exitCode == osqp::OsqpExitCode::kOptimal) {
            // get optimal solution
            const auto optimalObject = solver.objective_value();
            const auto optimalCoeffs = solver.primal_solution();

            //
            // format output:
            //
            for (int k = 0; k < K; ++k) {
                result.row(k) = optimalCoeffs.segment(k * N, N);
            }

            ROS_WARN("[Minimum Snap, Numeric]: Optimal objective is %.2f.", optimalObject);
        } else {
            // defaults to Eigen::MatrixXd::Zero():
            ROS_WARN("[Minimum Snap, Numeric]: Failed to find the optimal solution.");
        }
    } else {
        // defaults to Eigen::MatrixXd::Zero():
        ROS_WARN("[Minimum Snap, Numeric]: Failed to init OSQP solver.");
    }
    // done:
    return result;
}


Eigen::MatrixXd TrajectoryGeneratorWaypoint::PolyQPGeneration(
    const int d_order,           // the order of derivative
    const Eigen::MatrixXd &Path, // waypoints coordinates (3d)
    const Eigen::MatrixXd &Vel,  // boundary velocity
    const Eigen::MatrixXd &Acc,  // boundary acceleration
    const Eigen::VectorXd &Time) // time allocation in each segment
{
  // enforce initial and final velocity and accleration, for higher order
  // derivatives, just assume them be 0;
  int p_order = 2 * d_order - 1; // the order of polynomial
  int p_num1d = p_order + 1;     // the number of variables in each segment

  int m = Time.size();
  MatrixXd PolyCoeff(m, 3 * p_num1d);

  /**
   *
   * STEP 3.2:  generate a minimum-jerk piecewise monomial polynomial-based
   * trajectory
   *
   * **/
//    std::cout << "Vel Info : \n" << Vel << std::endl;
//    std::cout << "Acc Info : \n" << Acc << std::endl;
  // 使用osqp求解
//   ROS_INFO("d_order : %d,p_num1d : %d,m : %d, path(%d,%d),vel(%d,%d),acc(%d,%d)", d_order,p_num1d, m, Path.rows(),
//            Path.cols(), Vel.rows(), Vel.cols(), Acc.rows(), Acc.cols());
  
  for (int i = 0; i < Path.cols(); ++i) {
    PolyCoeff.block(0, i * p_num1d, m, p_num1d) =
        PolyQPGenerationNumeric(4, 3, Path.col(i),Vel.col(i), Acc.col(i), Time);
    }

  return PolyCoeff;
}

/**
 * @brief 获取目标函数的代价
 * 
 * @return double 
 */

double TrajectoryGeneratorWaypoint::getObjective() {
  _qp_cost = (_Px.transpose() * _Q * _Px + _Py.transpose() * _Q * _Py +
              _Pz.transpose() * _Q * _Pz)(0);
  return _qp_cost;
}

Vector3d TrajectoryGeneratorWaypoint::getPosPoly(MatrixXd polyCoeff, int k,
                                                 double t) {
  Vector3d ret;
  int _poly_num1D = (int)polyCoeff.cols() / 3;
  for (int dim = 0; dim < 3; dim++) {
    VectorXd coeff = (polyCoeff.row(k)).segment(dim * _poly_num1D, _poly_num1D);
    VectorXd time = VectorXd::Zero(_poly_num1D);

    for (int j = 0; j < _poly_num1D; j++)
      if (j == 0)
        time(j) = 1.0;
      else
        time(j) = pow(t, j);

    ret(dim) = coeff.dot(time);
    // cout << "dim:" << dim << " coeff:" << coeff << endl;
  }

  return ret;
}

Vector3d TrajectoryGeneratorWaypoint::getVelPoly(MatrixXd polyCoeff, int k,
                                                 double t) {
  Vector3d ret;
  int _poly_num1D = (int)polyCoeff.cols() / 3;
  for (int dim = 0; dim < 3; dim++) {
    VectorXd coeff = (polyCoeff.row(k)).segment(dim * _poly_num1D, _poly_num1D);
    VectorXd time = VectorXd::Zero(_poly_num1D);

    for (int j = 0; j < _poly_num1D; j++)
      if (j == 0)
        time(j) = 0.0;
      else
        time(j) = j * pow(t, j - 1);

    ret(dim) = coeff.dot(time);
  }

  return ret;
}

Vector3d TrajectoryGeneratorWaypoint::getAccPoly(MatrixXd polyCoeff, int k,
                                                 double t) {
  Vector3d ret;
  int _poly_num1D = (int)polyCoeff.cols() / 3;
  for (int dim = 0; dim < 3; dim++) {
    VectorXd coeff = (polyCoeff.row(k)).segment(dim * _poly_num1D, _poly_num1D);
    VectorXd time = VectorXd::Zero(_poly_num1D);

    for (int j = 0; j < _poly_num1D; j++)
      if (j == 0 || j == 1)
        time(j) = 0.0;
      else
        time(j) = j * (j - 1) * pow(t, j - 2);

    ret(dim) = coeff.dot(time);
  }

  return ret;
}