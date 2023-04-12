import numpy as np
from setuptools import setup, Extension
from Cython.Build import cythonize

extensions = []
extensions.append( Extension("CythonMaxCut",
        ["CythonMaxCut.pyx"],
        language="c++")
)

setup(
    ext_modules = cythonize(extensions,
        include_path = ["."] + [np.get_include()],
        annotate = True,
        language_level = "3"),
    zip_safe = False
)
