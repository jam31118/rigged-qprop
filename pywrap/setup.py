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

#libqprop = ('qprop', {
#    'sources': cc_source_filepaths,
#    'include_dirs':['../src/base'],
#    'extra_compile_args':[gxxoptions]}
#    )

ext_modules = [
    Extension( "_test",
        ["./src/pyc/_test.cc"] + cc_source_filepaths + ["../src/main/imag-prop.cc"],
        include_dirs=['./src/include','../src/base','../src/main'],
        extra_compile_args=gxx_args
        )
]

setup(
#    libraries = [libqprop],
    ext_modules=ext_modules
)

