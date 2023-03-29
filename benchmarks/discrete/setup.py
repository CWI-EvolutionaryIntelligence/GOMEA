import numpy as np
from setuptools import setup, Extension
from Cython.Build import cythonize
from pathlib import Path
import os
import gomea

extensions = []

libs_path = str(Path(gomea.__path__[0]).parent.absolute())
gomea_path = str(Path(gomea.__path__[0]).absolute())
libs_to_find = ['fitness', 'discrete', 'linkage', 'output']
libs_to_link = []
for root, dirs, files in os.walk(gomea_path):
    for lib_to_find in libs_to_find:
        for f in files:
            if f.startswith(lib_to_find) and f.endswith('.so'):
                libs_to_link.append("-l:"+str(f))
#print(libs_to_link)
# find names of library files within gomea_path to link with and add to 'libs_to_link'
extensions.append( Extension("CythonTrapFunction",
        ["CythonTrapFunction.pyx"],
        language="c++")
)
extensions.append( Extension("CythonTrapFunctionCDEF",
        ["CythonTrapFunction-cdef.pyx"],
        include_dirs=[libs_path],
        extra_compile_args=libs_to_link,
        extra_link_args=libs_to_link,
        library_dirs=[gomea_path],
        language="c++")
)

setup(
    ext_modules = cythonize(extensions,
        include_path = ["."] + [libs_path] + [np.get_include()],
        annotate = True,
        language_level = "3"),
    zip_safe = False
)
