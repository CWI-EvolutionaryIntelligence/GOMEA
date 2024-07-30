project('gomea', 'cpp', 'cython',
    default_options : ['c_std=c17', 'cpp_std=c++17'])

fs = import('fs')
py = import('python').find_installation(pure: false)

numpy_dep = dependency('numpy', required : false )
python_dep = py.dependency(embed: true)
cython = find_program('cython')

# Dependencies - as originally specified in setup.py.
# Note that we may want to move to using dependency() as above.
eigen_inc = include_directories('./gomea/lib/Eigen')
eigen_dep = declare_dependency(include_directories : eigen_inc)
cxxopts_inc = include_directories('./gomea/lib/cxxopts-3.1.1/include/')
cxxopts_dep = declare_dependency(include_directories : cxxopts_inc)

cython_args = [
    '--include-dir',
    meson.current_build_dir(),
    '-X always_allow_keywords=true'
]

# file listings.

# glob.glob("gomea/src/common/*.cpp") + glob.glob("gomea/src/utils/*.cpp")
common_src = files([
    'gomea/src/common/linkage_config.cpp',
    'gomea/src/common/linkage_model.cpp',
    'gomea/src/common/output_statistics.cpp',
    'gomea/src/common/partial_solution.cpp',
    'gomea/src/common/solution.cpp',

    'gomea/src/utils/embed.cpp',
    'gomea/src/utils/time.cpp',
    'gomea/src/utils/tools.cpp',
])
# glob.glob("gomea/src/fitness/*.cpp") + glob.glob("gomea/src/fitness/benchmarks-rv/*.cpp") + glob.glob("gomea/src/fitness/benchmarks-discrete/*.cpp")
fitness_src = files([
    'gomea/src/fitness/gbo_fitness.cpp',
    'gomea/src/fitness/py_bbo_fitness.cpp',
    'gomea/src/fitness/your_fitness_discrete.cpp',
    'gomea/src/fitness/your_fitness_realvalued.cpp',
    'gomea/src/fitness/bbo_fitness.cpp',
    'gomea/src/fitness/fitness.cpp',
    'gomea/src/fitness/py_gbo_fitness.cpp',
    'gomea/src/fitness/benchmarks-rv/SOREBChainStrong_t.cpp',
    'gomea/src/fitness/benchmarks-rv/circlesInASquareBBO_t.cpp',
    'gomea/src/fitness/benchmarks-rv/rosenbrockFunctionBBO_t.cpp',
    'gomea/src/fitness/benchmarks-rv/rosenbrockFunction_t.cpp',
    'gomea/src/fitness/benchmarks-rv/sphereFunction_t.cpp',
    'gomea/src/fitness/benchmarks-rv/SOREBChainStrongBBO_t.cpp',
    'gomea/src/fitness/benchmarks-rv/sphereFunctionBBO_t.cpp',
    'gomea/src/fitness/benchmarks-discrete/deceptiveTrap_t.cpp',
    'gomea/src/fitness/benchmarks-discrete/maxCut_t.cpp',
    'gomea/src/fitness/benchmarks-discrete/deceptiveTrapBBO_t.cpp',
    'gomea/src/fitness/benchmarks-discrete/maxCutBBO_t.cpp',
    'gomea/src/fitness/benchmarks-discrete/oneMax_t.cpp'
])

# glob.glob("gomea/src/discrete/*.cpp")
discrete_src = [
    'gomea/src/discrete/Config.cpp',
    'gomea/src/discrete/gomea.cpp',
    'gomea/src/discrete/gomeaIMS.cpp',
    'gomea/src/discrete/main.cpp',
    'gomea/src/discrete/Population.cpp',
    'gomea/src/discrete/utils.cpp'
]

# glob.glob("gomea/src/real_valued/*.cpp")
real_valued_src = [
    'gomea/src/real_valued/linkage_model.cpp',
    'gomea/src/real_valued/distribution.cpp',
    'gomea/src/real_valued/partial_solutionRV.cpp',
    'gomea/src/real_valued/solutionRV.cpp',
    'gomea/src/real_valued/rv-gomea.cpp',
    'gomea/src/real_valued/population.cpp',
    'gomea/src/real_valued/Config.cpp',
    'gomea/src/real_valued/tools.cpp'
]

py_dep = py.dependency()

gomea_incl = [include_directories('gomea', '.'), include_directories('.')]
# gomea_pxi = declare_dependency(sources: copy_ef, include_directories: gomea_incl)

# To handle cython header generation.
subdir('gomea')

# Common library symbols
gomea_common_lib = static_library(
    'gomealib_common',
    common_src,
    dependencies: [python_dep, numpy_dep],
    cython_args: cython_args,
)

# To share symbols for fitness across modules.
gomea_fitness_lib = static_library(
    'gomealib_fitness',
    ['gomea/fitness.pyx', gomea_fitness_cython] + fitness_src,
    override_options : ['cython_language=cpp'],
    link_with: [gomea_common_lib],
    dependencies: [python_dep, numpy_dep, cxxopts_dep, eigen_dep, gomea_fitness_cython_dep],
    cython_args: cython_args,
)

gomea_fitness_ex = py.extension_module(
    'fitness',
    ['gomea/fitness.pyx'],
    link_with: [gomea_common_lib, gomea_fitness_lib],
    override_options : ['cython_language=cpp'],
    dependencies: [numpy_dep, cxxopts_dep, eigen_dep, gomea_fitness_cython_dep],
    include_directories: [gomea_incl],
    cython_args: cython_args,
    install: true,
    subdir: 'gomea',
)

gomea_discrete_ex = py.extension_module(
    'discrete',
    ['gomea/discrete.pyx'] + discrete_src,
    link_with: [gomea_common_lib, gomea_fitness_lib],
    override_options : ['cython_language=cpp'],
    dependencies: [numpy_dep, cxxopts_dep, gomea_fitness_cython_dep],
    include_directories: [gomea_incl],
    cython_args: cython_args,
    install: true,
    subdir: 'gomea',
)

gomea_real_valued_ex = py.extension_module(
    'real_valued',
    ['gomea/real_valued.pyx', gomea_real_valued_cython] + real_valued_src,
    link_with: [gomea_common_lib, gomea_fitness_lib],
    override_options : ['cython_language=cpp'],
    dependencies: [numpy_dep, cxxopts_dep, eigen_dep, gomea_real_valued_cython_dep, gomea_fitness_cython_dep],
    include_directories: [gomea_incl],
    cython_args: cython_args,
    install: true,
    subdir: 'gomea',
)


gomea_linkage = py.extension_module(
    'linkage',
    ['gomea/linkage.pyx'],
    link_with: [gomea_common_lib],
    override_options : ['cython_language=cpp'],
    dependencies: [numpy_dep],
    include_directories: [gomea_incl],
    cython_args: cython_args,
    install: true,
    subdir: 'gomea',
)

gomea_output = py.extension_module(
    'output',
    ['gomea/output.pyx'],
    link_with: [gomea_common_lib],
    override_options : ['cython_language=cpp'],
    dependencies: [numpy_dep],
    include_directories: [gomea_incl],
    cython_args: cython_args,
    install: true,
    subdir: 'gomea',
)

py.install_sources([
        'gomea/__init__.py',
    ],
    pure: false,
    subdir: 'gomea',
)
