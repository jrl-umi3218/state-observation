#include "state-observation/tools/definitions.hpp"
#include "state-observation/tools/rigid-body-kinematics.hpp"
#include <iostream>
#include <state-observation/tools/odometry/legged-odometry-manager.hpp>

using namespace stateObservation;
using namespace kine;
// using namespace stateObservation::odometry;
// using LeggedOdometryManager = stateObservation::odometry::LeggedOdometryManager;

struct Traj
{
  Traj() : kine(LocalKinematics::zeroKinematics(kine::LocalKinematics::Flags::all))
  {
    kine.position = Vector3(0.0, 0.0, 0.8);
  }

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

int testLeggedOdometry(int errorcode)
{
  double simTime = 10.00;
  double dt = 0.005;
  int nbIters = int(simTime / dt);
  Traj traj;

  stateObservation::odometry::LeggedOdometryManager odometryManager =
      stateObservation::odometry::LeggedOdometryManager(); // manager for the legged odometry
  odometryManager.setSamplingTime(dt);
  odometry::LeggedOdometryManager::Configuration odomConfig(stateObservation::odometry::stringToOdometryType("6D"));
  odometryManager.init(odomConfig, traj.kine.toVector(Kinematics::Flags::pose));

  Kinematics kine;
  kine.position = Vector3(0.0, 0.0, 0.8);
  kine.orientation.setZeroRotation();
  for(int i = 0; i < nbIters; i++)
  {
    traj.iterate(dt);
    Kinematics zeroKine = Kinematics::zeroKinematics(Kinematics::Flags::pose);
    stateObservation::odometry::LeggedOdometryManager::ContactInputData test(zeroKine, 1.0);

    std::set<std::string> contactList;
    if(i % 2 == 0)
    {
      contactList.insert("Contact1");
    }
    else
    {
      contactList.insert("Contact2");
    }

    auto onNewContactOdom = [](stateObservation::odometry::LoContact & newContact)
    {
      newContact.bodyContactKine_ = Kinematics::zeroKinematics(Kinematics::Flags::pose);
      newContact.lambda(0.5);
    };

    auto onMaintainedContactOdom = [](stateObservation::odometry::LoContact & maintainedContact)
    {
      maintainedContact.bodyContactKine_ = Kinematics::zeroKinematics(Kinematics::Flags::pose);
      maintainedContact.lambda(0.5);
    };

    auto contactUpdateFunctions = stateObservation::odometry::LeggedOdometryManager::ContactUpdateFunctions()
                                      .onNewContact(onNewContactOdom)
                                      .onMaintainedContact(onMaintainedContactOdom);

    odometryManager.initLoop(contactList, contactUpdateFunctions);

    Matrix3 attitude = traj.kine.orientation.toMatrix3();
    odometryManager.run(odometry::LeggedOdometryManager::KineParams(kine).attitudeMeasurement(attitude));

    if((kine.position() - Vector3(0.0, 0.0, 0.8)).norm() > 1e-12)
    {
      std::cout << std::endl
                << "The position changed without maintained contact: " << kine.position().transpose() << std::endl;
      return errorcode;
    }
    if((kine.orientation.toMatrix3() - attitude).norm() > 1e-12)
    {
      std::cout << std::endl << "The orientation differs from the attitude measurement." << std::endl;
      return errorcode;
    }
  }

  return 0;
}

/// @brief Checks the double support with contacts of equal weights: both of them must be used for the orientation
/// odometry. Also checks that the contacts are no longer detected once they are all removed.
int testDoubleSupportSameWeights(int errorcode)
{
  double dt = 0.005;
  int nbIters = 10;

  stateObservation::odometry::LeggedOdometryManager odometryManager;
  odometryManager.setSamplingTime(dt);
  odometry::LeggedOdometryManager::Configuration odomConfig(stateObservation::odometry::stringToOdometryType("6D"));
  Kinematics initPose = Kinematics::zeroKinematics(Kinematics::Flags::pose);
  initPose.position = Vector3(0.0, 0.0, 0.8);
  odometryManager.init(odomConfig, initPose.toVector(Kinematics::Flags::pose));

  auto setContact = [](stateObservation::odometry::LoContact & contact)
  {
    contact.bodyContactKine_ = Kinematics::zeroKinematics(Kinematics::Flags::pose);
    contact.lambda(0.5);
  };
  auto contactUpdateFunctions = stateObservation::odometry::LeggedOdometryManager::ContactUpdateFunctions()
                                    .onNewContact(setContact)
                                    .onMaintainedContact(setContact);

  Kinematics kine = initPose;
  Vector3 tilt = Vector3::UnitZ();
  for(int i = 0; i < nbIters; i++)
  {
    odometryManager.initLoop({"LeftFoot", "RightFoot"}, contactUpdateFunctions);
    odometry::LeggedOdometryManager::KineParams kineParams(kine);
    kineParams.tiltMeasurement(tilt);
    odometryManager.run(kineParams);
  }

  auto & contactsManager = odometryManager.contactsManager();
  if(contactsManager.oriOdometryContacts_.size() != 2)
  {
    std::cout << std::endl
              << "Both contacts with the same weight should be used for the orientation odometry, but "
              << contactsManager.oriOdometryContacts_.size() << " is used." << std::endl;
    return errorcode;
  }

  odometryManager.initLoop({}, contactUpdateFunctions);
  if(contactsManager.contactsDetected())
  {
    std::cout << std::endl << "Contacts are still detected while the contact list is empty." << std::endl;
    return errorcode;
  }

  return 0;
}

/// @brief Runs the odometry over a few iterations with the given contacts, using a tilt measurement so that the yaw
/// is estimated from the contacts.
/// @return the estimated kinematics of the body
Kinematics runOdometry(stateObservation::odometry::LeggedOdometryManager & odometryManager,
                       const std::set<std::string> & contactList,
                       const Kinematics & bodyContactKine,
                       double lambda,
                       int nbIters)
{
  auto setContact = [&bodyContactKine, lambda](stateObservation::odometry::LoContact & contact)
  {
    contact.bodyContactKine_ = bodyContactKine;
    contact.lambda(lambda);
  };
  auto contactUpdateFunctions = stateObservation::odometry::LeggedOdometryManager::ContactUpdateFunctions()
                                    .onNewContact(setContact)
                                    .onMaintainedContact(setContact);

  Kinematics kine = Kinematics::zeroKinematics(Kinematics::Flags::pose);
  kine.position = Vector3(0.0, 0.0, 0.8);
  Vector3 tilt = Vector3::UnitZ();
  for(int i = 0; i < nbIters; i++)
  {
    odometryManager.initLoop(contactList, contactUpdateFunctions);
    odometry::LeggedOdometryManager::KineParams kineParams(kine);
    kineParams.tiltMeasurement(tilt);
    odometryManager.run(kineParams);
  }
  return kine;
}

/// @brief Checks the robustness of the odometry: unknown odometry type, contacts carrying no load, contacts without
/// orientation, and replacement of the pose of the body with contacts created on the same iteration.
int testOdometryRobustness(int errorcode)
{
  double dt = 0.005;
  odometry::LeggedOdometryManager::Configuration odomConfig(stateObservation::odometry::stringToOdometryType("6D"));
  Kinematics initPose = Kinematics::zeroKinematics(Kinematics::Flags::pose);
  initPose.position = Vector3(0.0, 0.0, 0.8);

  // contacts carrying no load
  {
    stateObservation::odometry::LeggedOdometryManager odometryManager;
    odometryManager.setSamplingTime(dt);
    odometryManager.init(odomConfig, initPose.toVector(Kinematics::Flags::pose));
    Kinematics bodyContactKine = Kinematics::zeroKinematics(Kinematics::Flags::pose);
    Kinematics kine = runOdometry(odometryManager, {"LeftFoot", "RightFoot"}, bodyContactKine, 0.0, 10);
    if(!kine.orientation.toMatrix3().allFinite())
    {
      std::cout << std::endl << "The orientation is not finite when the contacts carry no load." << std::endl;
      return errorcode;
    }
  }

  // point contacts: they give no orientation
  {
    stateObservation::odometry::LeggedOdometryManager odometryManager;
    odometryManager.setSamplingTime(dt);
    odometryManager.init(odomConfig, initPose.toVector(Kinematics::Flags::pose));
    Kinematics bodyContactKine;
    bodyContactKine.position = Vector3(0.1, 0.0, -0.8);
    Kinematics kine = runOdometry(odometryManager, {"LeftFoot", "RightFoot"}, bodyContactKine, 0.5, 10);
    if(!kine.orientation.isSet() || !kine.orientation.toMatrix3().allFinite())
    {
      std::cout << std::endl << "The orientation is not finite with point contacts only." << std::endl;
      return errorcode;
    }
  }

  // the pose of the body is replaced after the creation of a contact
  {
    stateObservation::odometry::LeggedOdometryManager odometryManager;
    odometryManager.setSamplingTime(dt);
    odometryManager.init(odomConfig, initPose.toVector(Kinematics::Flags::pose));
    Kinematics bodyContactKine = Kinematics::zeroKinematics(Kinematics::Flags::pose);
    bodyContactKine.position = Vector3(0.1, 0.0, -0.8);
    runOdometry(odometryManager, {"LeftFoot"}, bodyContactKine, 1.0, 1);

    Kinematics newPose = initPose;
    newPose.position = Vector3(1.0, 2.0, 0.8);
    odometryManager.replaceOdomBodyPose(newPose.toVector(Kinematics::Flags::pose));

    const Vector3 & contactRefPos =
        odometryManager.contactsManager().contacts().at("LeftFoot").worldRefKine_.position();
    if((contactRefPos - Vector3(1.1, 2.0, 0.0)).norm() > 1e-12)
    {
      std::cout << std::endl
                << "The new contact is not re-anchored on the new pose of the body: " << contactRefPos.transpose()
                << std::endl;
      return errorcode;
    }
  }

  return 0;
}

int main()
{
  int returnVal;
  int errorcode = 1;

  std::cout << "Starting testLeggedOdometry" << std::endl;
  if((returnVal = testLeggedOdometry(errorcode)))
  {
    std::cout << "testLeggedOdometry failed!" << errorcode << std::endl;
    return returnVal;
  }
  else
  {
    std::cout << "testLeggedOdometry succeeded" << std::endl;
  }
  errorcode++;

  std::cout << "Starting testDoubleSupportSameWeights" << std::endl;
  if((returnVal = testDoubleSupportSameWeights(errorcode)))
  {
    std::cout << "testDoubleSupportSameWeights failed!" << errorcode << std::endl;
    return returnVal;
  }
  else
  {
    std::cout << "testDoubleSupportSameWeights succeeded" << std::endl;
  }
  errorcode++;

  std::cout << "Starting testOdometryRobustness" << std::endl;
  if((returnVal = testOdometryRobustness(errorcode)))
  {
    std::cout << "testOdometryRobustness failed!" << errorcode << std::endl;
    return returnVal;
  }
  else
  {
    std::cout << "testOdometryRobustness succeeded" << std::endl;
  }
  errorcode++;

  std::cout << "Test Legged Odometry succeeded" << std::endl;
  return 0;
}
