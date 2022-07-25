from setuptools import Extension, setup
from Cython.Build import cythonize
import sys
import glob
import _version

if sys.version_info[0] == 2:
    raise Exception('Python 2.x is no longer supported')

with open("README.md", 'r') as f:
    long_description = f.read()

ext_discrete = Extension("DiscreteGOMEA",
        glob.glob("modules/discrete_gomea/src/cython/*.pyx") + glob.glob("modules/discrete_gomea/src/cpp/*.cpp"),
        include_dirs=["modules/discrete_gomea/include/"],
        language="c++",
        extra_compile_args=["-std=c++17"],
        extra_link_args=["-std=c++17"])

ext_real_valued = Extension("RealValuedGOMEA",
        glob.glob("modules/real_valued_gomea/src/cython/*.pyx") + glob.glob("modules/real_valued_gomea/src/cpp/*.cpp") +
        glob.glob("modules/utils/src/cpp/*.cpp"),
        include_dirs=["modules/real_valued_gomea/include/","modules/real_valued_gomea/src/cython/","modules/utils/include/"],
        language="c++",
        extra_compile_args=["-std=c++17"],
        extra_link_args=["-std=c++17"],
        libraries=["armadillo"])
        #libraries=["m", "armadillo", "blas", "lapack"]),
        #library_dirs=["../lib/boost_1_74_0/stage/lib/"],
    #extra_objects=["lib/boost_1_74_0/stage/lib/libboost_timer.a"])

ext_fitness = Extension("Fitness",
        glob.glob("modules/fitness/src/cython/*.pyx") + glob.glob("modules/fitness/src/cpp/*.cpp"),
        include_dirs=["modules/fitness/include/"],
        language="c++",
        extra_compile_args=["-std=c++17"],
        extra_link_args=["-std=c++17"])

#ext_utils = Extension("Utils",
#        glob.glob("modules/utils/src/cython/*.pyx") + glob.glob("modules/utils/src/cpp/*.cpp"),
#        include_dirs=["modules/utils/include/","modules/fitness/src/cython/"],
#        language="c++",
#        extra_compile_args=["-std=c++17"],
#        extra_link_args=["-std=c++17"])

extensions = [ext_fitness,ext_discrete,ext_real_valued]

setup(
    name="gomea",
    description='Library for the use of various variants of the Gene-pool Optimal Mixing Evolutionary Algorith (GOMEA).',
    author='Anton Bouter',
    author_email='anton.bouter@cwi.nl',
    url='',
    version=_version.__version__,
    long_description=long_description,
    long_description_content_type='text/markdown',
    ext_modules=cythonize(extensions,include_path=glob.glob("modules/*/src/cython/"),gdb_debug=False,language_level = "3"),
    zip_safe=False
)

