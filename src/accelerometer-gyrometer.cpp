#include <state-observation/sensors-simulation/accelerometer-gyrometer.hpp>

namespace stateObservation
{
AccelerometerGyrometer::AccelerometerGyrometer(bool matrixMode, bool withAcceleroBias, bool withGyroBias)
: r_(Matrix3::Zero()), acc_(Vector3::Zero()), omega_(Vector3::Zero()), output_(Vector::Zero(measurementSize_, 1))
{
#ifdef STATEOBSERVATION_VERBOUS_CONSTRUCTORS
  std::cout << std::endl << "AccelerometerGyrometer Constructor" << std::endl;
#endif // STATEOBSERVATION_VERBOUS_CONSTRUCTOR
  withGyroBias_ = withGyroBias;
  withAcceleroBias_ = withAcceleroBias;
  matrixMode_ = matrixMode;
  updateStateSize_();
}

Index AccelerometerGyrometer::getStateSize_() const
{
  return currentStateSize_;
}

Index AccelerometerGyrometer::getMeasurementSize_() const
{
  return measurementSize_;
}

Vector AccelerometerGyrometer::computeNoiselessMeasurement_()
{

  int biasIndex;
  if(!matrixMode_)
  {
    Quaternion q(state_.head<4>());

    r_ = q.toRotationMatrix();
    acc_ = state_.segment<3>(4);
    omega_ = state_.segment<3>(7);
    biasIndex = 10;
  }
  else
  {
    r_ = Eigen::Map<Matrix3>(&state_[0]);
    acc_ = state_.segment<3>(9);
    omega_ = state_.segment<3>(12);
    biasIndex = 15;
  }

  output_.head<3>() = accelerationMeasure(acc_, r_);
  output_.tail<3>() = rotationVelocityMeasure(omega_, r_);

  /// Add a bias
  if(withGyroBias_)
  {
    output_.head<3>().noalias() += state_.segment<3>(biasIndex);
    biasIndex += 3;
  }

  if(withAcceleroBias_)
  {
    output_.tail<3>().noalias() += state_.segment<3>(biasIndex);
  }

  return output_;
}

void AccelerometerGyrometer::updateStateSize_()
{
  Index accBiasSize, gyroBiasSize;
  if(withAcceleroBias_)
  {
    accBiasSize = 3;
  }
  else
  {
    accBiasSize = 0;
  }

  if(withGyroBias_)
  {
    gyroBiasSize = 3;
  }
  else
  {
    gyroBiasSize = 0;
  }

  if(!matrixMode_)
  {
    currentStateSize_ = stateSize_ + accBiasSize + gyroBiasSize;
  }
  else
  {
    currentStateSize_ = stateSizeMatrix_ + accBiasSize + gyroBiasSize;
  }
}

} // namespace stateObservation
