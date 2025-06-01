/**
 * \file      tilt-estimator.hpp
 * \author    Rafael Cisneros, Mehdi Benallegue
 * \date       2018
 * \brief      Defines the class for the tilt estimator.
 *
 * \details
 *
 *
 */

#ifndef TILTESTIMATORHPP
#define TILTESTIMATORHPP

#include <state-observation/api.h>
#include <state-observation/observer/zero-delay-observer.hpp>

namespace stateObservation
{

/**
 * \class  TiltEstimator
 * \brief
 *         Description is pending
 *
 *         use getEstimatedState to obtain the state vector
 *         the tilt R.transpose()*e_z is constituted
 *         with the last three components of the state vector.

           This estimator is based on the work detailed in:
            Benallegue, Mehdi, et al. "Lyapunov-stable orientation estimator for humanoid robots." IEEE Robotics and
            Automation Letters 5.4 (2020): 6371-637
 *
 */
class STATE_OBSERVATION_DLLAPI TiltEstimator : public ZeroDelayObserver
{
public:
  /// The constructor
  ///  \li alpha : parameter related to the convergence of the linear velocity
  ///              of the IMU expressed in the control frame
  ///  \li beta  : parameter related to the fast convergence of the tilt
  ///  \li gamma : parameter related to the orthogonality
  ///  \li dt    : sampling time
  TiltEstimator(double alpha, double beta, double gamma, double dt);

  /// @brief initializes the state vector.
  /// @param x1 The initial local linear velocity of the IMU.
  /// @param x2_p The initial value of the intermediate estimate of the IMU's tilt.
  /// @param x2 The initial tilt of the IMU.
  void initEstimator(Vector3 x1 = Vector3::Zero(), Vector3 x2_prime = Vector3::UnitZ(), Vector3 x2 = Vector3::UnitZ());

  /// set the gain of x1_hat variable
  void setAlpha(const double alpha)
  {
    alpha_ = alpha;
  }
  double getAlpha() const
  {
    return alpha_;
  }

  /// set the gain of x2prime_hat variable
  void setBeta(const double beta)
  {
    beta_ = beta;
  }
  double getBeta() const
  {
    return beta_;
  }

  /// set the gain of x2_hat variable
  void setGamma(const double gamma)
  {
    gamma_ = gamma;
  }
  double getGamma() const
  {
    return gamma_;
  }

  /// set the sampling time of the measurements
  void setSamplingTime(const double dt)
  {
    dt_ = dt;
  }
  double getSamplingTime() const
  {
    return dt_;
  }

  /// sets ths measurement (accelero and gyro stacked in one vector)
  void setMeasurement(const Vector3 & yv_k, const Vector3 & ya_k, const Vector3 & yg_k, TimeIndex k);

  using ZeroDelayObserver::setMeasurement;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
  /// The parameters of the estimator
  double alpha_, beta_, gamma_;

  /// Sampling time
  double dt_;

  /// variables used for the computation
  Vector3 x1_;
  Vector3 x1_hat_;
  Vector3 x2_hat_prime_;
  Vector3 x2_hat_;
  Vector3 dx1_hat;

  /// The tilt estimator loop
  StateVector oneStepEstimation_();
};

} // namespace stateObservation

#endif // TILTESTIMATORHPP
