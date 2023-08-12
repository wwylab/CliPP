from setuptools import setup, Extension
import sys
import os

#if not (sys.version_info[0] == 3 and sys.version_info[1] >= 5 and sys.version_info[2] >= 1):
version_morph = sys.version_info[0]*10000+sys.version_info[1]*100+sys.version_info[2]
version_base = 30501
if not (version_morph >= version_base):
    sys.stderr.write("Error message: CliP can only run with python >=3.5.1\n")
    sys.exit(-1)

if sys.platform.startswith('darwin'):
    os.environ['CC'] = "clang"
    os.environ['CXX'] = "clang++"
    ext_modules=[Extension('CliP', ['./src/kernel.cpp'], include_dirs=['eigen-3.4.0'], extra_compile_args = ['-O3'], extra_link_args = ['-O3']),]
else:
    ext_modules=[Extension('CliP', ['./src/kernel.cpp'], include_dirs=['eigen-3.4.0'], extra_compile_args = ['-O3', '-fopenmp'], extra_link_args = ['-O3', '-fopenmp']),]

    
setup(
    name="CliP",
    ext_modules=ext_modules)


