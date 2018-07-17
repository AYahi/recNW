from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
  name = 'RecNW_aligner',
  ext_modules=[
    Extension(name="recnw",
              sources=["recnw_cython.pyx", "recnw.cpp"], # Note, you can link against a c++ library instead of including the source
              language="c++")
              # extra_link_args=['-framework', 'OpenGL', '-framework', 'GLUT']),
    ],
  cmdclass = {'build_ext': build_ext},

)
