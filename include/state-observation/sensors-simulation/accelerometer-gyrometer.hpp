/**
 * \file      accelerometer-gyrometer.hpp
 * \author    Mehdi Benallegue
 * \date       2013
 * \brief     Implements the accelerometer-gyrometer inertial measuremen
 *
 *
 */

#ifndef SIMULATIONACCELEROMETERGYROMETERSENSORHPP
#define SIMULATIONACCELEROMETERGYROMETERSENSORHPP

#include <boost/assert.hpp>
#include <Eigen/Core>

#include <state-observation/api.h>
#include <state-observation/sensors-simulation/algebraic-sensor.hpp>
#include <state-observation/sensors-simulation/algorithm/linear-acceleration.hpp>
#include <state-observation/sensors-simulation/algorithm/rotation-velocity.hpp>

namespace stateObservation
{
/**
 * \class  AccelerometerGyrometer
 * \brief  Implements the accelerometer-gyrometer measurements
 *
 *
 *
 * \details
 *
 */

class STATE_OBSERVATION_DLLAPI AccelerometerGyrometer : public AlgebraicSensor,
                                                        protected algorithm::LinearAcceleration,
                                                        protected algorithm::RotationVelocity
{
public:
  AccelerometerGyrometer(bool matrixMode = false, bool withAcceleroBias = false, bool withGyroBias = false);

  /// Virtual destructor
  virtual ~AccelerometerGyrometer() {}

  void setMatrixMode(bool matrixMode)
  {
    matrixMode_ = matrixMode;
    updateStateSize_();
  }

  void setWithGyroBias(bool withGyroBias)
  {
    withGyroBias_ = withGyroBias;
    updateStateSize_();
  }

  void setWithAcceleroBias(bool withAcceleroBias)
  {
    withAcceleroBias_ = withAcceleroBias;
    updateStateSize_();
  }

protected:
  /// Gets the state vector Size
  virtual Index getStateSize_() const;

  /// Gets the measurements vector size
  virtual Index getMeasurementSize_() const;

  virtual Vector computeNoiselessMeasurement_();

  void updateStateSize_();

  Matrix3 r_;
  Vector3 acc_;
  Vector3 omega_;
  Vector output_;

  bool withGyroBias_;
  bool withAcceleroBias_;

  bool matrixMode_;

  static const Index stateSize_ = 10;
  static const Index stateSizeMatrix_ = 15;
  static const Index measurementSize_ = 6;

  Index currentStateSize_;
};

} // namespace stateObservation

#endif // SIMULATIONACCELEROMETERGYROMETERSENSORHPP
