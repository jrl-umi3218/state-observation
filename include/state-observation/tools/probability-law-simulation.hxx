template<typename BiasType, typename StdType>
Eigen::Matrix<double, BiasType::RowsAtCompileTime, BiasType::ColsAtCompileTime>
    ProbabilityLawSimulation::getGaussianMatrix(BiasType bias, StdType std, Index rows, Index cols)
{
  typedef Eigen::Matrix<double, BiasType::RowsAtCompileTime, BiasType::ColsAtCompileTime> ReturnType;

  static_assert(isEigen<StdType>::value && isEigen<BiasType>::value,
                "Standard deviation and bias need to be Eigen matrices");

  std::normal_distribution<double> g(0, 1);

  if(rows == -1)
  {
    rows = bias.rows();
  }
  if(cols == -1)
  {
    rows = bias.cols();
  }
  return bias + std * FixOrDynMatrixTools<ReturnType>::nullaryExp([&]() { return getGaussianScalar(); }, rows, cols);
}

template<typename ReturnType>
typename MatrixType<ReturnType>::type ProbabilityLawSimulation::getUniformMatrix(Index rows,
                                                                                 Index cols,
                                                                                 double min,
                                                                                 double max)
{
  return FixOrDynMatrixTools<ReturnType>::nullaryExp([&]() { return getUniformScalar(min, max); }, rows, cols);
}