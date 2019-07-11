from distutils.core import setup, Extension

module1 = Extension('ccmap',
                    libraries = ['m'],
                    include_dirs = ['./include'],
                    sources = ['ccmapmodule.c', './src/mesh.c', './src/transform_mesh.c', './src/decoygen.c'],
                    extra_compile_args=['-D', 'AS_PYTHON_EXTENSION', '-std=c99', '-pedantic'])
                    #extra_compile_args=['-D', 'DEBUG', '-D', 'AS_PYTHON_EXTENSION', '-std=c99', '-pedantic'])

setup (name = 'ccmapModule',
       version = '1.0',
       description = 'This is the C implementation of the mesh based contact map',
       ext_modules = [module1])
