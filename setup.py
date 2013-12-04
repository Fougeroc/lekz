from distutils.core import setup
from distutils.extension import Extension

try:
    import Cython
    USE_CYTHON = True
    print "using cython"
except ImportError:
    USE_CYTHON = False
    print "cython not present, use the default C file"

extension = ".pyx" if USE_CYTHON else ".c"

sourcefiles = [
# cython files
"src/lekz" + extension,
# pure C file
"src/generalized_permutation.c",
"src/lin_alg.c",
"src/quad_cyclic_cover.c",
"src/random.c"]

extensions = [Extension("lekz", sourcefiles)]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

setup(
    name="lekz",
    version = "0.alpha1",
    author = "V. Delecroix",
    author_email = "vincent.delecroix@labri.fr",
    description="Computing Lyapunov exponents of the Kontsevich-Zorich cocycle",
    long_description=
"""LEKZ: Lyapunov exponents of the Kontzevich-Zorich cocycles

Each affine SL(2,R)-ergodic measure on the moduli space of quadratic
differentials on Riemman surfaces give rise to Lyapunov exponents. This package
allows to compute approximations of them for various affine measures:
- the Masur-Veech measure
- coverings construction (which contain the so called pillowcase covers and
  origamis)""",
    ext_modules = extensions)
