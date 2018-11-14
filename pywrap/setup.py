from distutils.core import setup, Extension

from os import listdir
from os.path import splitext, join
base_dirpath = '../src/base'
ll = listdir(base_dirpath)
cc_source_files = [fn for fn in ll if splitext(fn)[-1] == '.cc']
cc_source_files.remove('kbhit.cc')
cc_source_filepaths = [join(base_dirpath, fn) for fn in cc_source_files]

gxx_args = ["-g","-std=c++0x","-O8","-Wno-deprecated","-march=native",
        "-funroll-all-loops","-finline-functions","-fexpensive-optimizations",
        "-ffast-math","-Xlinker","-defsym","-Xlinker","MAIN__=main","-I."
        ]

from os import environ
QPROP_DEP_DIR = environ.get("QPROP_DEP_DIR")
from os.path import isdir
if not isdir(QPROP_DEP_DIR):
    raise IOError("QPROP_DEP_DIR({}) could not be found".format(QPROP_DEP_DIR))
GSL_HOME = join(QPROP_DEP_DIR, "gsl")

#libqprop = ('qprop', {
#    'sources': cc_source_filepaths,
#    'include_dirs':['../src/base'],
#    'extra_compile_args':[gxxoptions]}
#    )

sources_filepaths = ["./c/_test.cc"] \
        + cc_source_filepaths \
        + ["../src/main/potentials.cc"] \
        + ["../src/main/imag-prop.cc", "../src/main/real-prop.cc", "../src/main/eval-tsurff.cc"]

ext_modules = [
    Extension( "_test",
        sources_filepaths,
        include_dirs=['../src/base','../src/main', join(GSL_HOME, "include")],
        extra_compile_args=gxx_args,
        library_dirs=[join(GSL_HOME,"lib")],
        libraries=["gsl","gslcblas"]
        )
]

setup(
#    libraries = [libqprop],
    ext_modules=ext_modules
)

