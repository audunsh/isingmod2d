from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize('isingmod.pyx'),  # accepts a glob pattern
)