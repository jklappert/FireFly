
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

test('bunch', test_bunch)
test('normal', test_normal)
test('safe', test_safe)
test('save read', test_save_read)
test('save write', test_save_write)
