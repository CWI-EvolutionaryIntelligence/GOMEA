VERSION_MAJOR = 0
VERSION_MINOR = 0
VERSION_PATCH = 1
VERSION_STRING = '%s.%s.%s' % (VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH)
__version__ = VERSION_STRING

from setuptools import Extension, setup, find_namespace_packages
from Cython.Build import cythonize
import sys
import glob
import numpy as np

if sys.version_info[0] == 2:
    raise Exception('Python 2.x is no longer supported')

with open("README.md", 'r') as f:
    long_description = f.read()

common_src = glob.glob("gomea/src/common/*.cpp") + glob.glob("gomea/src/utils/*.cpp")
fitness_src = glob.glob("gomea/src/fitness/*.cpp") + glob.glob("gomea/src/fitness/benchmarks-rv/*.cpp")

extensions = []

extensions.append( Extension("gomea.discrete",
        ["gomea/discrete.pyx"] + glob.glob("gomea/src/discrete/*.cpp") + common_src + fitness_src,
        include_dirs=["."],
        language="c++",
        extra_compile_args=["-std=c++17"],
        extra_link_args=["-std=c++17"])
)

extensions.append( Extension("gomea.real_valued",
        ["gomea/real_valued.pyx"] + glob.glob("gomea/src/real_valued/*.cpp") + common_src + fitness_src,
        include_dirs=["."],
        language="c++",
        extra_compile_args=["-std=c++17"],
        extra_link_args=["-std=c++17"],
        libraries=["armadillo"],
        library_dirs=[],
        extra_objects=[])
)

extensions.append( Extension("gomea.fitness",
        ["gomea/fitness.pyx"] + fitness_src + common_src,
        include_dirs=["."],
        language="c++",
        extra_compile_args=["-std=c++17"],
        extra_link_args=["-std=c++17"])
)

setup(
    name = "gomea",
    description = 'Library for the use of various variants of the Gene-pool Optimal Mixing Evolutionary Algorith (GOMEA).',
    author = 'Anton Bouter',
    author_email = 'Anton.Bouter@cwi.nl',
    url = '',
    version = __version__,
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    packages=["gomea"],
    ext_modules = cythonize(extensions,
        include_path = ["."] + [np.get_include()],
        gdb_debug = False,
        language_level = "3"),
    zip_safe = False
)

