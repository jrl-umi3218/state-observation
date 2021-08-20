#include <state-observation/tools/probability-law-simulation.hpp>

namespace stateObservation
{
namespace tools
{
const int ProbabilityLawSimulation::defaultSeed_ = 1985;
std::random_device ProbabilityLawSimulation::rd_;
std::mt19937 ProbabilityLawSimulation::gen_(ProbabilityLawSimulation::defaultSeed_);

double ProbabilityLawSimulation::getGaussianScalar(double std, double bias)
{
  std::normal_distribution<> g(bias, std);
  return g(gen_);
}

double ProbabilityLawSimulation::getUniformScalar(double min, double max)
{
  std::uniform_real_distribution<double> g(min, max);
  return g(gen_);
}

void ProbabilityLawSimulation::setSeed(unsigned int seed)
{
  gen_.seed(seed);
}

void ProbabilityLawSimulation::setRandomSeed()
{
  gen_.seed(ProbabilityLawSimulation::rd_());
}
} // namespace tools

} // namespace stateObservation
