project('discrete_frechet',
        'cpp',
        default_options: ['cpp_std=c++20',
                          'warning_level=3',
                          'werror=true'],
        version: '0.0.1')

cpp = meson.get_compiler('cpp')

src_inc_dir = include_directories('src')

test_exec = executable('frechet_test',
                       ['src/test.cpp'],
                       include_directories: src_inc_dir)
