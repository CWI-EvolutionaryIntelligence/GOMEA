# Manually cythonize so that we get the header file `fitness.h`
gomea_fitness_cython = custom_target(
    'fitness.pyx',
    output : ['fitness.h'],
    input : ['fitness.pyx'],
    command : [cython, '-3', '--cplus', '@INPUT@', '-o', '@OUTPUT@'],
)
gomea_fitness_cython_dep = declare_dependency(
    sources: [gomea_fitness_cython],
)

# Manually cythonize so that we get the header file `real_valued.h`
gomea_real_valued_cython = custom_target(
    'real_valued.pyx',
    output : ['real_valued.h'],
    input : 'real_valued.pyx',
    command : [cython, '-3', '--cplus', '@INPUT@', '-o', '@OUTPUT@'],
)
gomea_real_valued_cython_dep = declare_dependency(
    sources: [gomea_real_valued_cython],
)

copy_ef = fs.copyfile('EmbeddedFitness.pxi')