#include "state-observation/observer/tilt-estimator-humanoid.hpp"
#include "state-observation/observer/tilt-estimator.hpp"
#include "state-observation/tools/definitions.hpp"
#include "state-observation/tools/rigid-body-kinematics.hpp"
#include <iostream>

using namespace stateObservation;
using namespace kine;

struct Traj
{
  Traj() : kine(LocalKinematics::zeroKinematics(kine::LocalKinematics::Flags::all)) {}

  void iterate(double dt, bool inMotion = true)
  {
    Vector3 linJerk = tools::ProbabilityLawSimulation::getUniformMatrix(3, 1, -0.2, 0.2) * 1;
    Vector3 angJerk = tools::ProbabilityLawSimulation::getUniformMatrix(3, 1, -0.2, 0.2) * 1;

    if(inMotion)
    {
      kine.linAcc() += linJerk * dt;
      kine.angAcc() += angJerk * dt;
    }
    else
    {
      double K = 20;
      kine.linAcc() = -K * kine.linVel();
      kine.angAcc() = -K * kine.angVel();
    }

    kine.integrate(dt);
  }
  Vector3 getX2()
  {
    return kine.orientation.toMatrix3().transpose() * Vector3::UnitZ();
  }
  Vector getYa()
  {
    return kine.linAcc() + cst::gravityConstant * getX2();
  }
  LocalKinematics & operator()()
  {
    return kine;
  }
  LocalKinematics kine;
};

int testTiltObserver(int errorcode)
{
  double simTime = 10.00;
  double dt = 0.005;
  int nbIters = int(simTime / dt);

  double err;
  TiltEstimator tilt(5, 2, 2, dt);
  Traj traj;

  Vector3 x2_init = traj.getX2();
  tilt.initEstimator(traj().linVel(), x2_init, x2_init);

  for(int i = 0; i < nbIters; i++)
  {
    traj.iterate(dt);

    tilt.setMeasurement(traj().linVel(), traj.getYa(), traj().angVel(), i + 1);
    tilt.getEstimatedState(i + 1);

    err = (tilt.getCurrentEstimatedState().segment<3>(0) - traj().linVel()).squaredNorm();
    if(err > 1e-4)
    {
      std::cout << std::endl << "The local velocity estimate is incorrect." << std::endl;
      std::cout << std::endl << "Estimated: " << tilt.getCurrentEstimatedState().segment<3>(0).transpose() << std::endl;
      std::cout << std::endl << "Simulated: " << traj().linVel().transpose() << std::endl;
      return errorcode;
    }

    err = (tilt.getCurrentEstimatedState().segment<3>(6) - traj.getX2()).squaredNorm();
    if(err > 1e-4)
    {
      std::cout << std::endl << "The tilt estimate is incorrect." << std::endl;
      std::cout << std::endl << "Estimated: " << tilt.getCurrentEstimatedState().segment<3>(6).transpose() << std::endl;
      std::cout << std::endl << "Simulated: " << traj.getX2().transpose() << std::endl;

      return errorcode;
    }
  }

  return 0;
}

/// @brief Checks that the estimator gives the same result when it has to catch up several iterations at once (all the
/// measurements are given before calling getEstimatedState) as when it is updated at every iteration.
int testMultiStepCatchUp(int errorcode)
{
  double simTime = 1.00;
  double dt = 0.005;
  int nbIters = int(simTime / dt);
  int catchUpSteps = 4;

  TiltEstimator tiltStepByStep(5, 2, 2, dt);
  TiltEstimator tiltCatchUp(5, 2, 2, dt);
  Traj traj;

  Vector3 x2_init = traj.getX2();
  tiltStepByStep.initEstimator(traj().linVel(), x2_init, x2_init);
  tiltCatchUp.initEstimator(traj().linVel(), x2_init, x2_init);

  for(int i = 0; i < nbIters; i++)
  {
    traj.iterate(dt);

    tiltStepByStep.setMeasurement(traj().linVel(), traj.getYa(), traj().angVel(), i + 1);
    tiltCatchUp.setMeasurement(traj().linVel(), traj.getYa(), traj().angVel(), i + 1);

    tiltStepByStep.getEstimatedState(i + 1);

    if((i + 1) % catchUpSteps == 0)
    {
      tiltCatchUp.getEstimatedState(i + 1);

      double err = (tiltCatchUp.getCurrentEstimatedState() - tiltStepByStep.getCurrentEstimatedState()).squaredNorm();
      if(tiltCatchUp.getCurrentTime() != i + 1 || err > 1e-20)
      {
        std::cout << std::endl
                  << "The estimate obtained by catching up " << catchUpSteps
                  << " iterations differs from the step-by-step one." << std::endl;
        std::cout << "Step by step: " << tiltStepByStep.getCurrentEstimatedState().transpose() << std::endl;
        std::cout << "Catch up:     " << tiltCatchUp.getCurrentEstimatedState().transpose() << std::endl;
        return errorcode;
      }
    }
  }

  return 0;
}

/// @brief Checks the measurement of TiltEstimatorHumanoid computed from the kinematics of the control frame in the IMU
/// frame: it must be the local linear velocity of the IMU in the world when the control frame is fixed in the world.
int testHumanoidMeasurement(int errorcode)
{
  TiltEstimatorHumanoid tilt(5, 2, 2, 0.005);
  Vector3 x1 = Vector3::Zero();
  Vector3 x2 = Vector3::UnitZ();
  Vector3 x2_prime = Vector3::UnitZ();
  tilt.initEstimator(x1, x2_prime, x2);

  for(int i = 0; i < 10; i++)
  {
    // kinematics of the IMU in the world
    Matrix3 R = kine::Orientation::randomRotation().toMatrix3();
    Vector3 worldImuPos = tools::ProbabilityLawSimulation::getUniformMatrix<Vector3>();
    Vector3 worldImuLinVel = tools::ProbabilityLawSimulation::getUniformMatrix<Vector3>();
    Vector3 yg = tools::ProbabilityLawSimulation::getUniformMatrix<Vector3>();
    // the control frame is fixed in the world
    Vector3 worldControlPos = tools::ProbabilityLawSimulation::getUniformMatrix<Vector3>();

    // kinematics of the control frame in the IMU frame
    Vector3 imuControlPos = R.transpose() * (worldControlPos - worldImuPos);
    Vector3 imuControlLinVel = -yg.cross(imuControlPos) - R.transpose() * worldImuLinVel;

    tilt.setMeasurement(imuControlPos, imuControlLinVel, Vector3::Zero(), yg, i + 1);

    double err = (tilt.getMeasurement(i + 1).head<3>() - R.transpose() * worldImuLinVel).norm();
    if(err > 1e-12)
    {
      std::cout << std::endl << "The local linear velocity of the IMU is incorrect, error: " << err << std::endl;
      return errorcode;
    }
  }

  return 0;
}

int main()
{
  int returnVal;
  int errorcode = 1;

  std::cout << "Starting testTiltObserver" << std::endl;
  if((returnVal = testTiltObserver(errorcode)))
  {
    std::cout << "testTiltObserver failed!" << errorcode << std::endl;
    return returnVal;
  }
  else
  {
    std::cout << "testTiltObserver succeeded" << std::endl;
  }
  errorcode++;

  std::cout << "Starting testMultiStepCatchUp" << std::endl;
  if((returnVal = testMultiStepCatchUp(errorcode)))
  {
    std::cout << "testMultiStepCatchUp failed!" << errorcode << std::endl;
    return returnVal;
  }
  else
  {
    std::cout << "testMultiStepCatchUp succeeded" << std::endl;
  }
  errorcode++;

  std::cout << "Starting testHumanoidMeasurement" << std::endl;
  if((returnVal = testHumanoidMeasurement(errorcode)))
  {
    std::cout << "testHumanoidMeasurement failed!" << errorcode << std::endl;
    return returnVal;
  }
  else
  {
    std::cout << "testHumanoidMeasurement succeeded" << std::endl;
  }
  errorcode++;

  std::cout << "Test Tilt Observer succeeded" << std::endl;
  return 0;
}
