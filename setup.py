from distutils.core import setup, Extension
import os
import sysconfig

extra_compile_args = sysconfig.get_config_var('CFLAGS').split()
extra_compile_args += ['-D', 'AS_PYTHON_EXTENSION', '-std=c99', '-pedantic']
module1 = Extension('ccmap',
                    libraries = ['m'],
                    include_dirs = ['./include'],#, '../modules/include/python3.6m/'],
                    sources = ['ccmapmodule.c', 'src/decoygen.c', 'src/encode.c', 'src/molecular_object.c', 'src/cell_crawler.c', 'src/mesh.c', 'src/transform_mesh.c', 'src/mesh_default.c', 'src/pdb_coordinates.c'],
                    extra_compile_args=extra_compile_args
		    #extra_compile_args=['-D', 'DEBUG', '-D', 'AS_PYTHON_EXTENSION', '-std=c99', '-pedantic'])
		)

setup (name = 'ccmapModule',
       version = '1.0',
       description = 'This is the C implementation of the mesh based contact map',
       ext_modules = [module1])


