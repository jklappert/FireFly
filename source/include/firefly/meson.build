
firefly_version = meson.project_version().split('.')

version_data = configuration_data()
version_data.set('FireFly_VERSION_MAJOR', firefly_version[0])
version_data.set('FireFly_VERSION_MINOR', firefly_version[1])
version_data.set('FireFly_VERSION_RELEASE', firefly_version[2])
version_target = configure_file(
  output : 'version.hpp',
  configuration : version_data)

config_data = configuration_data()
config_data.set('FLINT', cpp_def_flint)
config_data.set('DEFAULT', cpp_def_default)
config_data.set('CUSTOM', false)
config_data.set('WITH_MPI', cpp_def_mpi)
config_target = configure_file(
  output : 'config.hpp',
  configuration : config_data)

headers = [
  'AmplitudeParser.hpp',
  'BaseReconst.hpp',
  'BlackBoxBase.hpp',
  'DenseSolver.hpp',
  'FFInt.hpp',
  'FFIntVec.hpp',
  'FFThieleInterpolator.hpp',
  'gzstream.hpp',
  'HornerGenerator.hpp',
  'Logger.hpp',
  'Monomial.hpp',
  'ParserUtils.hpp',
  'Poly.hpp',
  'Polynomial.hpp',
  'PolynomialFF.hpp',
  'PolyReconst.hpp',
  'RationalFunction.hpp',
  'RationalFunctionFF.hpp',
  'RationalNumber.hpp',
  'RatReconst.hpp',
  'ReconstHelper.hpp',
  'Reconstructor.hpp',
  'ShuntingYardParser.hpp',
  'ThreadPool.hpp',
  'tinydir.h',
  'UintHasher.hpp',
  'utils.hpp']

headers_mpi = ['MPIWorker.hpp']

install_headers(headers, subdir : 'firefly')
install_headers(version_target, subdir : 'firefly')

if cpp_def_mpi
  install_headers(headers_mpi, subdir : 'firefly')
endif
