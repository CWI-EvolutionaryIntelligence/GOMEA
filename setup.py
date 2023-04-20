VERSION_MAJOR = 1
VERSION_MINOR = 0
VERSION_PATCH = 0
VERSION_STRING = '%s.%s.%s' % (VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH)
__version__ = VERSION_STRING

from setuptools import Extension, setup, find_namespace_packages
from Cython.Build import cythonize
import sys
import glob
import numpy as np
import platform

if sys.version_info[0] == 2:
        raise Exception('Python 2.x is no longer supported')

with open("README.md", 'r') as f:
        long_description = f.read()

debug_mode = False
if '--debug' in sys.argv:
        debug_mode = True

common_src = glob.glob("gomea/src/common/*.cpp") + glob.glob("gomea/src/utils/*.cpp")
fitness_src = glob.glob("gomea/src/fitness/*.cpp") + glob.glob("gomea/src/fitness/benchmarks-rv/*.cpp") + glob.glob("gomea/src/fitness/benchmarks-discrete/*.cpp")

compile_args = ["-std=c++17"]
link_args = ["-std=c++17"]
if debug_mode:
        compile_args.extend(['-UNDEBUG','-g'])
        link_args.extend(['-UNDEBUG','-g'])
else:
        compile_args.extend(['-O3','-g0'])
        link_args.extend(['-O3','-g0'])
if platform.system() == "Darwin":
        compile_args.extend(["-stdlib=libc++","-mmacosx-version-min=10.15"])
        link_args.extend(["-stdlib=libc++","-mmacosx-version-min=10.15"])

if platform.system() == "Windows":
        compile_args = ["/std:c++17"]
        link_args = ["/std:c++17"]
        if debug_mode:
                compile_args.extend(['/UNDEBUG','/Zi','/g0'])
                link_args.extend(['/UNDEBUG','/Zi','/g0'])
        else:
                compile_args.extend(['/O2'])
                link_args.extend(['/O2'])

extensions = []

extensions.append( Extension("gomea.discrete",
        ["gomea/discrete.pyx"] + glob.glob("gomea/src/discrete/*.cpp") + common_src + fitness_src,
        include_dirs = [".", "lib/cxxopts-3.1.1/include/"] + [np.get_include()],
        language="c++",
        extra_compile_args=compile_args,
        extra_link_args=link_args)
)

extensions.append( Extension("gomea.real_valued",
        ["gomea/real_valued.pyx"] + glob.glob("gomea/src/real_valued/*.cpp") + common_src + fitness_src,
        include_dirs = ["."] + ["lib/Eigen"] + [np.get_include()],
        language="c++",
        extra_compile_args=compile_args,
        extra_link_args=link_args,
        library_dirs=[],
        extra_objects=[])
)

extensions.append( Extension("gomea.fitness",
        ["gomea/fitness.pyx"] + fitness_src + common_src,
        include_dirs = ["."] + [np.get_include()],
        language="c++",
        extra_compile_args=compile_args,
        extra_link_args=link_args)
)

extensions.append( Extension("gomea.linkage",
        ["gomea/linkage.pyx"] + common_src,
        include_dirs = ["."] + [np.get_include()],
        language="c++",
        extra_compile_args=compile_args,
        extra_link_args=link_args)
)

extensions.append( Extension("gomea.output",
        ["gomea/output.pyx"] + common_src,
        include_dirs = ["."] + [np.get_include()],
        language="c++",
        extra_compile_args=compile_args,
        extra_link_args=link_args)
)

setup(
    name = "gomea",
    description = 'Library for the use of various variants of the Gene-pool Optimal Mixing Evolutionary Algorith (GOMEA).',
    author = 'Anton Bouter',
    author_email = 'Anton.Bouter@cwi.nl',
    url = 'https://github.com/abouter/gomea',
    version = __version__,
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    packages=["gomea"],
    include_dirs=["gomea"],
    ext_modules = cythonize(extensions,
        include_path = ["."] + [np.get_include()],
        gdb_debug = debug_mode,
        language_level = "3"),
    install_requires=["numpy>=1.23.0","tqdm"],
    include_package_data=True,
    package_data = {
        'gomea': ['*.pxd', '*.hpp', '*.h']
    },
    zip_safe = False
)

