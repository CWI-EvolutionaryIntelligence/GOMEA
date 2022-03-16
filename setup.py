#from setuptools import setup

#from Cython.Build import cythonize

#setup(ext_modules=cythonize("src/pyx/*.pyx"),include_path=["inc/hpp/"])
from setuptools import Extension, setup
from Cython.Build import cythonize
import sys
import _version

if sys.version_info[0] == 2:
    raise Exception('Python 2.x is no longer supported')

with open("README.md", 'r') as f:
    long_description = f.read()

extensions = [
    Extension("pyrvgomea", ["src/rvgomea/cython/pyrvgomea.pyx"],
        include_dirs=["inc/rvgomea/hpp/","src/rvgomea/cpp/"])
        #include_dirs=["lib/boost_1_74_0/","inc/hpp/","src/cpp/"],
        #library_dirs=["../lib/boost_1_74_0/stage/lib/"],
        #libraries=["boost_timer"])
        # Libraries need to be installed locally    
    #extra_objects=["lib/boost_1_74_0/stage/lib/libboost_timer.a"])
    # Everything but the above is included here.
    #Extension("*", ["*.pyx"],
    #    include_dirs=[""],
]
setup(
    name="gomea",
    description='Library for the use of various variants of the Model-Based Evolutionary Algorith GOMEA.',
    author='Anton Bouter',
    author_email='anton.bouter@cwi.nl',
    url='',
    version=_version.__version__,
    long_description=long_description,
    long_description_content_type='text/markdown',
    ext_modules=cythonize(extensions),
    zip_safe=False
)

