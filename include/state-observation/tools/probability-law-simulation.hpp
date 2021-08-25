/**
 * \file      probability-law-simulation.hpp
 * \author    Mehdi Benallegue
 * \date       2013
 * \brief
 *
 *
 *
 */

#ifndef SENSORSSIMULATIONPROBABILITYLAWSIMULATIONHPP
#define SENSORSSIMULATIONPROBABILITYLAWSIMULATIONHPP

#include <random>

#include <state-observation/api.h>
#include <state-observation/tools/definitions.hpp>

namespace stateObservation
{
namespace tools
{
class STATE_OBSERVATION_DLLAPI ProbabilityLawSimulation
{
public:
  /// gets a scalar Gaussian random variable
  /// having a given bias and standard deviation(std)
  /// default is the cetered unit Gaussian
  static double getGaussianScalar(double std = 1, double bias = 0);

  /// gets vector Gaussian random variable
  /// having a given bias and standard deviation(std)
  template<typename BiasType,
           typename StdType = Eigen::Matrix<double, BiasType::RowsAtCompileTime, BiasType::RowsAtCompileTime>>
  static Eigen::Matrix<double, BiasType::RowsAtCompileTime, BiasType::ColsAtCompileTime> getGaussianMatrix(
      BiasType bias = Eigen::Matrix<double, BiasType::RowsAtCompileTime, BiasType::ColsAtCompileTime>::Zero(),
      StdType std = Eigen::Matrix<double, BiasType::RowsAtCompileTime, BiasType::RowsAtCompileTime>::Identity(),
      Index rows = Index(EigenType<BiasType>::type::RowsAtCompileTime),
      Index cols = Index(EigenType<BiasType>::type::ColsAtCompileTime));

  /// @brief Get a scalar following a uniform distribution between min and max
  ///
  /// @param min the minimal value of the variable
  /// @param max the maximal value of the variable
  /// @return double The simulated random variable
  static double getUniformScalar(double min = -1., double max = 1.);

  /// @brief Get a Matrix following a uniform distribution between min and max
  ///
  /// @param min the minimal value of the variable
  /// @param max the maximal value of the variable
  template<typename ReturnType = Matrix>
  static typename MatrixType<ReturnType>::type getUniformMatrix(Index rows = ReturnType::RowsAtCompileTime,
                                                                Index cols = ReturnType::ColsAtCompileTime,
                                                                double min = -1.,
                                                                double max = 1.);

  /// @brief sets the seed to the generator
  static void setSeed(unsigned int seed);

  /// @brief sets the seed to the generator
  static void setRandomSeed();

protected:
  static const int defaultSeed_;
  /// @brief Private constructor for preventing instanciation of Probability Law Simulation
  ProbabilityLawSimulation() {}
  static std::random_device rd_;
  static std::mt19937 gen_;
};
#include <state-observation/tools/probability-law-simulation.hxx>
} // namespace tools
} // namespace stateObservation
#endif // SENSORSSIMULATIONPROBABILITYLAWSIMULATIONHPP
