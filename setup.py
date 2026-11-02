"""Build hook for the hopacity wheel.

Compiles libhop.so with standalone/Makefile (gfortran by default; COMPILER=ifort or ifx
in the environment selects those) and places it in the package next to the header and
HOP's default data files.  The header is a link to standalone/libhop.h in the source
tree and a file in the wheel.  The library links the compiler's Fortran and
OpenMP runtimes dynamically, so a wheel for other machines has to bundle them
(delocate on macOS, auditwheel on Linux).
"""
import shutil
import subprocess
from pathlib import Path

from setuptools import Distribution, setup
from setuptools.command.bdist_wheel import bdist_wheel
from setuptools.command.build_py import build_py

HERE = Path(__file__).resolve().parent


class BuildWithLibrary(build_py):
    def run(self):
        super().run()
        pkg = Path(self.build_lib) / 'hopacity'
        build = Path(self.get_finalized_command('build').build_temp).resolve() / 'libhop'
        lib = build / 'libhop.so'
        subprocess.run(['make', '-C', str(HERE / 'standalone'), f'BUILD={build}', str(lib)],
                       check=True)
        shutil.copyfile(lib, pkg / 'libhop.so')
        shutil.copyfile(HERE / 'hopacity' / 'libhop.h', pkg / 'libhop.h')
        src, dst = HERE / 'hopacity' / 'data', pkg / 'data'
        shutil.copytree(src, dst, dirs_exist_ok=True)


class BinaryDistribution(Distribution):
    """Carries a compiled library, so the wheel is platform-specific."""

    def has_ext_modules(self):
        return True


class PlatformWheel(bdist_wheel):
    """A wheel for this platform but any Python 3: the library is loaded through cffi's
    ABI mode, so it does not depend on the interpreter's C API."""

    def get_tag(self):
        return 'py3', 'none', super().get_tag()[2]


setup(distclass=BinaryDistribution,
      cmdclass={'build_py': BuildWithLibrary, 'bdist_wheel': PlatformWheel})
