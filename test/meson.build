
test_bunch = executable(
  'test_bunch',
  ['test_bunch.cpp'],
  include_directories : include_dirs,
  dependencies : deps,
  cpp_args : rootdir_cppdef,
  link_with : firefly,
  install : false)

test_normal = executable(
  'test_normal',
  ['test_normal.cpp'],
  include_directories : include_dirs,
  dependencies : deps,
  cpp_args : rootdir_cppdef,
  link_with : firefly,
  install : false)

test_safe = executable(
  'test_safe',
  ['test_safe.cpp'],
  include_directories : include_dirs,
  dependencies : deps,
  cpp_args : rootdir_cppdef,
  link_with : firefly,
  install : false)

test_save_read = executable(
  'test_save_read',
  ['test_save_read.cpp'],
  include_directories : include_dirs,
  dependencies : deps,
  cpp_args : rootdir_cppdef,
  link_with : firefly,
  install : false)

test_save_write = executable(
  'test_save_write',
  ['test_save_write.cpp'],
  include_directories : include_dirs,
  dependencies : deps,
  cpp_args : rootdir_cppdef,
  link_with : firefly,
  install : false)

test('bunch', test_bunch,
  env : 'MALLOC_PERTURB_=0',
  is_parallel : false)

test('normal', test_normal,
  env : 'MALLOC_PERTURB_=0',
  is_parallel : false)

test('safe', test_safe,
  env : 'MALLOC_PERTURB_=0',
  is_parallel : false)

test('save write', test_save_write,
  env : 'MALLOC_PERTURB_=0',
  is_parallel : false)

test('save read', test_save_read,
  env : 'MALLOC_PERTURB_=0',
  is_parallel : false)
