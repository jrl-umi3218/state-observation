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

#ifndef TILTESTIMATORHUMANOIDHPP
#define TILTESTIMATORHUMANOIDHPP

#include <state-observation/api.h>
#include <state-observation/observer/tilt-estimator.hpp>
#include <state-observation/observer/zero-delay-observer.hpp>

namespace stateObservation
{

/**
 * \class  TiltEstimatorHumanoid
 * \brief
 *         Description is pending
 *
 *         use getEstimatedState to obtain the state vector
 *         the tilt R.transpose()*e_z is constituted
 *         with the last three components of the state vector.
 *
 */
class STATE_OBSERVATION_DLLAPI TiltEstimatorHumanoid : public TiltEstimator
{
public:
  /// The constructor
  ///  \li alpha : parameter related to the convergence of the linear velocity
  ///              of the IMU expressed in the control frame
  ///  \li beta  : parameter related to the fast convergence of the tilt
  ///  \li gamma : parameter related to the orthogonality
  TiltEstimatorHumanoid(double alpha, double beta, double gamma);

  /// sets the position of the IMU sensor in the control frame
  void setSensorPositionInC(const Vector3 & p)
  {
    p_S_C_ = p;
  }

  Vector3 getSensorPositionInC()
  {
    return p_S_C_;
  }

  /// sets the oriantation of the IMU sensor in the control frame
  void setSensorOrientationInC(const Matrix3 & R)
  {
    R_S_C_ = R;
  }
  Matrix3 getSensorOrientationInC()
  {
    return R_S_C_;
  }

  Vector3 getVirtualLocalVelocityMeasurement()
  {
    return x1_;
  }

  /// sets teh linear velocity of the IMU sensor in the control frame
  void setSensorLinearVelocityInC(const Vector3 & v)
  {
    v_S_C_ = v;
  }

  Vector3 getSensorLinearVelocityInC()
  {
    return v_S_C_;
  }

  /// sets the angular velocity of the IMU sensor in the control frame
  void setSensorAngularVelocityInC(const Vector3 & w)
  {
    w_S_C_ = w;
  }
  Vector3 getSensorAngularVelocityInC()
  {
    return w_S_C_;
  }

  /// sets the velocity of the control origin in the world frame
  /// this velocity has to be expressed in the control frame.
  void setControlOriginVelocityInW(const Vector3 & v)
  {
    v_C_ = v;
  }
  Vector3 getControlOriginVelocityInW()
  {
    return v_C_;
  }

  /// @brief informs the estimator that x1hat (the estimate of the local linear velocity of the IMU in the world) needs
  /// to be reset.
  /// @copydetails checkResetX1hat()
  void resetImuLocVelHat();

/// prevent c++ overloaded virtual function warning
#if defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Woverloaded-virtual"
#else
#  if defined(__GNUC__)
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Woverloaded-virtual"
#  endif
#endif

  // we also want to use the function setMeasurement from the TiltEstimator class, that is hidden by the following
  using TiltEstimator::setMeasurement;
  /// sets ths measurement (accelero and gyro stacked in one vector)
  void setMeasurement(const Vector3 & ya_k, const Vector3 & yg_k, TimeIndex k);

#if defined(__clang__)
#  pragma clang diagnostic pop
#else
#  if defined(__GNUC__)
#    pragma GCC diagnostic pop
#  endif
#endif

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

public:
protected:
  /// Position of the IMU in the control frame
  Vector3 p_S_C_;

  /// Orientation of the IMU in the control frame
  Matrix3 R_S_C_;

  /// Linear velocity of the IMU in the control frame
  Vector3 v_S_C_;

  /// Angular velocity of the IMU in the control frame
  Vector3 w_S_C_;

  /// Linear velocity of the control frame
  Vector3 v_C_;

  /// indicates that x1hat needs to be reset
  bool resetX1hat_ = false;

  /// @brief Checks if x1hat needs to be reset and if yes, resets it with the current value of x1.
  /// @details x1hat needs to be reset when the mode of computation of the anchor frame changed, to avoid
  /// discontinuities.
  void checkResetX1hat();
};

} // namespace stateObservation

#endif // TILTESTIMATORHUMANOIDHPP
