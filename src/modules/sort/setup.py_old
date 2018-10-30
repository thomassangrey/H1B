from distutils.core import setup, Extension
from distutils.sysconfig import get_python_inc

c_ext = Extension("sort_map", ["sort_map.c"])

setup(
    ext_modules=[c_ext],
    include_dirs=[get_python_inc(), '../utilities'] 
)