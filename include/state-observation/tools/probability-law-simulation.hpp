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
  template<typename ReturnType = Matrix, typename StdType, typename BiasType>
  static typename MatrixType<ReturnType>::type getGaussianMatrix(StdType std,
                                                                 BiasType bias,
                                                                 Index rows = BiasType::RowsAtCompileTime,
                                                                 Index cols = BiasType::ColsAtCompileTime);

  /// @brief Get a scalar following a uniform distribution between min and max
  ///
  /// @param min the minimal value of the variable
  /// @param max the maximal value of the variable
  /// @return double The simulated random variable
  static double getUniformScalar(double min = 0., double max = 1.);

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
