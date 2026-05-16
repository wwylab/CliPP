from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import os
from pathlib import Path

#if not (sys.version_info[0] == 3 and sys.version_info[1] >= 5 and sys.version_info[2] >= 1):
version_morph = sys.version_info[0]*10000+sys.version_info[1]*100+sys.version_info[2]
version_base = 30501
if not (version_morph >= version_base):
    sys.stderr.write("Error message: CliPP can only run with python >=3.5.1\n")
    sys.exit(-1)

def _nvidia_package_paths(module_name):
    import importlib.util

    spec = importlib.util.find_spec(module_name)
    if spec is None or not spec.submodule_search_locations:
        return None

    root = Path(next(iter(spec.submodule_search_locations)))
    include_dir = root / "include"
    lib_dir = root / "lib"
    if not include_dir.is_dir() or not lib_dir.is_dir():
        return None
    return include_dir, lib_dir


def _find_libcuda_dir():
    import ctypes.util

    found = ctypes.util.find_library("cuda")
    if found:
        found_path = Path(found)
        if found_path.is_absolute() and found_path.parent.is_dir():
            return found_path.parent

    for candidate in (
        "/usr/lib/x86_64-linux-gnu",
        "/usr/lib64",
        "/usr/local/cuda/lib64",
        "/usr/local/cuda/lib",
    ):
        candidate_path = Path(candidate)
        if (candidate_path / "libcuda.so").exists() or (candidate_path / "libcuda.so.1").exists():
            return candidate_path

    return None


def _env_flag(name):
    value = os.environ.get(name)
    if value is None or value == "":
        return None

    normalized = value.lower()
    if normalized in {"1", "true", "yes", "on"}:
        return True
    if normalized in {"0", "false", "no", "off"}:
        return False

    sys.stderr.write("Error message: %s must be one of 1/0, true/false, yes/no, or on/off.\n" % name)
    sys.exit(-1)


def _detect_cuda_build():
    libcuda_dir = _find_libcuda_dir()
    if libcuda_dir is None:
        return {
            "available": False,
            "runtime_paths": None,
            "nvrtc_paths": None,
            "libcuda_dir": None,
            "missing": ["libcuda"],
        }

    cuda_runtime_paths = _nvidia_package_paths("nvidia.cuda_runtime")
    cuda_nvrtc_paths = _nvidia_package_paths("nvidia.cuda_nvrtc")

    missing = []
    if cuda_runtime_paths is None:
        missing.append("nvidia-cuda-runtime-cu12")
    if cuda_nvrtc_paths is None:
        missing.append("nvidia-cuda-nvrtc-cu12")
    return {
        "available": len(missing) == 0,
        "runtime_paths": cuda_runtime_paths,
        "nvrtc_paths": cuda_nvrtc_paths,
        "libcuda_dir": libcuda_dir,
        "missing": missing,
    }


requested_cuda = _env_flag("CLIPP_USE_CUDA")
cuda_info = None
use_cuda = False

if requested_cuda is True:
    cuda_info = _detect_cuda_build()
    if not cuda_info["available"]:
        sys.stderr.write(
            "Error message: CUDA build requested but CUDA build dependencies were not found. "
            "Missing: %s.\n" % ", ".join(cuda_info["missing"])
        )
        sys.exit(-1)
    use_cuda = True
elif requested_cuda is None:
    cuda_info = _detect_cuda_build()
    use_cuda = cuda_info["available"]

include_dirs = ['eigen-3.4.0']
extra_compile_args = ['-O3', '-std=c++17']
extra_link_args = ['-O3']
define_macros = []
library_dirs = []
libraries = []
runtime_library_dirs = []
sources = ['./src/kernel_cpu.cpp']

if use_cuda:
    if cuda_info is None:
        cuda_info = _detect_cuda_build()

    cuda_runtime_include, cuda_runtime_lib = cuda_info["runtime_paths"]
    cuda_nvrtc_include, cuda_nvrtc_lib = cuda_info["nvrtc_paths"]
    libcuda_dir = cuda_info["libcuda_dir"]

    include_dirs.extend([str(cuda_runtime_include), str(cuda_nvrtc_include)])
    import importlib.util

    triton_spec = importlib.util.find_spec("triton")
    if triton_spec is not None and triton_spec.origin:
        triton_cuda_include = Path(triton_spec.origin).parent / "backends" / "nvidia" / "include"
        if triton_cuda_include.is_dir():
            include_dirs.append(str(triton_cuda_include))
    define_macros.append(("USE_CUDA", "1"))
    libraries.append("cuda")
    library_dirs.append(str(libcuda_dir))
    runtime_library_dirs.extend([str(cuda_runtime_lib), str(cuda_nvrtc_lib)])
    extra_link_args.extend([
        str(cuda_runtime_lib / "libcudart.so.12"),
        str(cuda_nvrtc_lib / "libnvrtc.so.12"),
    ])
    sources.append('./src/kernel_cuda.cpp')

if sys.platform.startswith('darwin'):
    os.environ['CC'] = "clang"
    os.environ['CXX'] = "clang++"
    ext_modules=[Extension('CliPP', sources, include_dirs=include_dirs, define_macros=define_macros, extra_compile_args=extra_compile_args, extra_link_args=extra_link_args),]
else:
    ext_modules=[Extension('CliPP', sources, include_dirs=include_dirs, define_macros=define_macros, library_dirs=library_dirs, libraries=libraries, runtime_library_dirs=runtime_library_dirs, extra_compile_args=extra_compile_args + ['-fopenmp'], extra_link_args=extra_link_args + ['-fopenmp']),]


class BuildExt(build_ext):
    def finalize_options(self):
        super().finalize_options()
        self.force = True

    
setup(
    name="CliPP",
    ext_modules=ext_modules,
    cmdclass={"build_ext": BuildExt})