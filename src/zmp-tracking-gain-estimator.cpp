#include <state-observation/dynamics-estimators/zmp-tracking-gain-estimator.hpp>

namespace stateObservation
{
/// definitions of static variables
constexpr double ZmpTrackingGainEstimator::defaultGainDriftSecond;
constexpr double ZmpTrackingGainEstimator::defaultZmpProcessErrorStd;
constexpr double ZmpTrackingGainEstimator::defaultZmpMeasurementErrorStd;
constexpr double ZmpTrackingGainEstimator::defaultGainMinimum;
constexpr double ZmpTrackingGainEstimator::defaultZmpUncertainty;
constexpr double ZmpTrackingGainEstimator::defaultGainUncertainty;

ZmpTrackingGainEstimator::ZmpTrackingGainEstimator(double dt,
                                                   const Vector2 & zmpMeasureErrorStd,
                                                   const Vector3 & gainDriftPerSecond,
                                                   const Vector2 & zmpProcessErrorStd,
                                                   double minimumGain,
                                                   const Vector2 & initZMP,
                                                   const Vector3 & initGain,
                                                   const Vector2 & initZMPUncertainty,
                                                   const Vector3 & initGainUncertainty)
: dt_(dt), minimumGain_(minimumGain), zmpMeasureErrorstd_(zmpMeasureErrorStd),
  gainDriftPerSecondStd_(gainDriftPerSecond), zmpProcessErrorStd_(zmpProcessErrorStd), filter_(5, 2, 0)
{
  A_.setIdentity();
  C_.topLeftCorner<2, 2>().setIdentity();
  C_.bottomRightCorner<2, 3>().setZero();
  filter_.setC(C_);
  Q_.setZero();
  resetWithMeasurements(initZMP, initGain, 0., initZMPUncertainty, initGainUncertainty);
}

void ZmpTrackingGainEstimator::resetWithMeasurements(const Vector2 & initZMP,
                                                     const Vector3 & initGain,
                                                     const Matrix2 & yaw,
                                                     const Vector2 & initZMPUncertainty,
                                                     const Vector3 & initGainUncertainty)
{
  filter_.reset();
  updateQ_();
  updateR_();

  Vector5 x;
  x << initZMP, initGain;
  filter_.setState(x, 0);

  yaw_ = yaw;

  Matrix5 P;
  // clang-format off
  P << Vec2ToSqDiag_(initZMPUncertainty),  Matrix23::Zero(),
       Matrix23::Zero().transpose(),       Vec3ToSqDiag_(initGainUncertainty);
  // clang-format on
  filter_.setStateCovariance(P);
}

void ZmpTrackingGainEstimator::setSamplingTime(double dt)
{
  dt_ = dt;
}

void ZmpTrackingGainEstimator::setGain(const Vector3 & gain)
{
  Vector5 x = filter_.getCurrentEstimatedState();
  // update the gain part of the state
  x.tail<3>() = gain;
  filter_.setCurrentState(x);
}

void ZmpTrackingGainEstimator::setGain(const Vector3 & gain, const Vector3 & uncertainty)
{
  setGain(gain);
  Matrix5 P = filter_.getStateCovariance();
  /// resetting the non diagonal parts
  P.topRightCorner<2, 3>().setZero();
  P.bottomLeftCorner<3, 2>().setZero();
  /// updating the gain part
  P.bottomRightCorner<3, 3>() = Vec3ToSqDiag_(uncertainty);
  filter_.setStateCovariance(P);
}

void ZmpTrackingGainEstimator::setGainDriftPerSecond(const Vector3 & gaindrift)
{
  Q_.bottomRightCorner<3, 3>() = Vec3ToSqDiag_(gaindrift);
  filter_.setProcessCovariance(Q_);
}

void ZmpTrackingGainEstimator::setMinimumGain(const double & minGain)
{
  minimumGain_ = minGain;
}

void ZmpTrackingGainEstimator::setZMPProcesError(const Vector2 & /*zmpprocesserrorstd*/) {}

/// @brief Set the Zmp Measurement Error Stamdard deviation
///
void ZmpTrackingGainEstimator::setZmpMeasureErrorStd(const Vector2 &) {}

Matrix2 ZmpTrackingGainEstimator::getEstimatedLocalGain() const
{
  Vector5 x = filter_.getCurrentEstimatedState();
  Matrix2 K;

  // clang-format off
  K<<x(2),x(4),
     x(4),x(3);
  // clang-format on
  return K;
}

/// @brief Set the Inputs of the estimator.
///
/// @param zmperr  measurement of the DCM in the world frame
/// @param zmp  mesaurement of the ZMP in the world frame
/// @param R    the 2x2 Matrix'representing the yaw angle i.e. bias_global == R * bias*local
void ZmpTrackingGainEstimator::setInputs(const Vector2 & /*zmpErr*/, const Vector2 & /*zmp*/, const Matrix2 & /*R*/) {}

/// @brief Runs the estimation. Needs to be called every timestep
///
/// @return Vector2
void ZmpTrackingGainEstimator::update() {}

} // namespace stateObservation
