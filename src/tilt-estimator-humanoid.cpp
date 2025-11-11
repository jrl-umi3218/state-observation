#include <state-observation/observer/tilt-estimator-humanoid.hpp>

namespace stateObservation
{
TiltEstimatorHumanoid::TiltEstimatorHumanoid() : TiltEstimator() {}

TiltEstimatorHumanoid::TiltEstimatorHumanoid(double alpha, double beta, double gamma, double dt)
: TiltEstimator(alpha, beta, gamma, dt)
{
}

TiltEstimatorHumanoid::TiltEstimatorHumanoid(double alpha, double beta, double gamma, int n, int m, double dt)
: TiltEstimator(alpha, beta, gamma, n, m, dt)
{
}

void TiltEstimatorHumanoid::setMeasurement(const Vector3 & imuControlPos,
                                           const Vector3 & imuControlLinVel,
                                           const Vector3 & ya_k,
                                           const Vector3 & yg_k,
                                           TimeIndex k)
{
  auto start = std::chrono::high_resolution_clock::now();
  x1_ = -yg_k.cross(imuControlPos) - imuControlLinVel;

  auto end = std::chrono::high_resolution_clock::now();

  iterTime_ += std::chrono::duration<double, std::micro>(end - start).count();

  TiltEstimator::setMeasurement(x1_, ya_k, yg_k, k);
}

void TiltEstimatorHumanoid::resetImuLocVelHat()
{
  x_().segment<3>(0) = x1_;
}

} // namespace stateObservation
