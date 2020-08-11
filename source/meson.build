
src = [
  'AmplitudeParser.cpp',
  'BaseReconst.cpp',
  'FFInt.cpp',
  'FFIntVec.cpp',
  'FFThieleInterpolator.cpp',
  'gzstream.cpp',
  'HornerGenerator.cpp',
  'Monomial.cpp',
  'ParserUtils.cpp',
  'Poly.cpp',
  'Polynomial.cpp',
  'PolynomialFF.cpp',
  'PolyReconst.cpp',
  'RationalFunction.cpp',
  'RationalFunctionFF.cpp',
  'RationalNumber.cpp',
  'RatReconst.cpp',
  'ReconstHelper.cpp',
  'ShuntingYardParser.cpp',
  'utils.cpp']

if get_option('static')
  firefly = static_library(
    'firefly',
    src,
    include_directories : include_dirs,
    dependencies : deps,
    install_rpath : rpath,
    install : true)
else
  firefly = library(
    'firefly',
    src,
    include_directories : include_dirs,
    dependencies : deps,
    install_rpath : rpath,
    install : true)
endif

pkg = import('pkgconfig')
pkg.generate(
  firefly,
  filebase : 'firefly',
  name : 'FireFly',
  description :
    'Library for the reconstruction of rational functions over finite fields.',
  requires : pkgconfig_deps,
  libraries : other_deps)

# Dependency object to be used if FireFly is built as a subproject:
firefly_dep = declare_dependency(
  include_directories : include_dirs,
  dependencies : pkgconfig_deps + other_deps,
  link_with : firefly)

if meson.version().version_compare('>=0.54.0')
  meson.override_dependency('firefly', firefly_dep)
endif
