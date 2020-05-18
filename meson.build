
project(
  'firefly',
  ['cpp'],
  version : '2.0.1',
  license : 'GPL3',
  default_options : [
    'cpp_std=c++14',
    'buildtype=release',
    'warning_level=3',
    'default_library=both'])

if get_option('arch') != ''
  add_project_arguments('-march=' + get_option('arch'), language : ['cpp'])
endif

include_dirs = [include_directories('source/include')]
# pkgconfig_deps will be listed in the 'Requires' field of the pkgfconfig file.
pkgconfig_deps = [dependency('zlib')]
# other_deps will be listed in the 'Libs' field
# (with flags in the 'Cflags' field).
other_deps = [dependency('threads')]
cpp_defs = []
rootdir_cppdef = '-DFIREFLY_ROOT_DIR="' + meson.current_source_dir() + '"'
rpath = ''

# Meson versions before 0.49 don't provide access to the option wrap_mode.
wrap_mode = 'default'
if meson.version().version_compare('>=0.49.0')
  wrap_mode = get_option('wrap_mode')
endif

cpp_compiler = meson.get_compiler('cpp')

# On some systems, gmp doesn't provide a pkg-config file.
gmp_dep = dependency('gmp', required : false)
if gmp_dep.found()
  pkgconfig_deps += gmp_dep
else
  gmp_dep = cpp_compiler.find_library('gmp', required : true)
  other_deps += gmp_dep
endif

if get_option('flint')
  cpp_defs += ['-DFLINT']
  if wrap_mode == 'forcefallback'
    flint_found = false
  else
    # FLint doesn't provide a pkg-config file.
    # This is also the reason why we can't use
    # dependency('flint', fallback : ['flint', 'flint_dep']))
    flint_dep = cpp_compiler.find_library('flint', required : false)
    flint_found = flint_dep.found()
  endif
  if not flint_found
    if wrap_mode == 'nofallback'
      error('FLINT subproject build is disabled.')
    endif
    flint_proj = subproject('flint2')
    flint_dep = flint_proj.get_variable('flint_dep')
    rpath = join_paths(get_option('prefix'),get_option('libdir'))
  endif
  other_deps += [flint_dep]
else
  warning(
    'Using default implementation of modular arithmetic. ' +
    'This can be much slower than FLINT.')
  cpp_defs += ['-DDEFAULT']
endif

if get_option('mpi')
  cpp_defs += ['-DWITH_MPI']
  pkgconfig_deps += [dependency('mpi', language : 'cpp')]
  example_src = ['example_mpi.cpp']
else
  example_src = ['example.cpp']
endif

if get_option('jemalloc')
  pkgconfig_deps += dependency('jemalloc')
endif

add_project_arguments(cpp_defs, language : 'cpp')
deps = pkgconfig_deps + other_deps

test_link_args = []
if cpp_compiler.has_link_argument('-Wl,--disable-new-dtags')
  test_link_args = ['-Wl,--disable-new-dtags']
endif

subdir('source')
subdir('source/include/firefly')
subdir('test')

executable(
  'example',
  example_src,
  include_directories : include_dirs,
  link_with : firefly,
  dependencies : deps,
  cpp_args : rootdir_cppdef,
  link_args : test_link_args,
  install : false)

executable(
  'ff_insert',
  'ff_insert.cpp',
  include_directories : include_dirs,
  link_with : firefly,
  dependencies : deps,
  install_rpath : rpath,
  install : true)
