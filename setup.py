from distutils.core import setup, Extension
import os
import sysconfig

extra_compile_args = sysconfig.get_config_var('CFLAGS').split()
extra_compile_args += ['-D', 'AS_PYTHON_EXTENSION', '-std=c99', '-pedantic']
<<<<<<< HEAD
extra_link_args = ['-L/home/julia/IBCP/modules/lib/python3.6/site-packages/']
module1 = Extension('ccmap',
                    libraries = ['m'],
                    include_dirs = ['./include', '../modules/include/python3.6m/'],
                    sources = ['ccmapmodule.c', './src/mesh.c', './src/transform_mesh.c', './src/decoygen.c','./src/encode.c'],
=======
print(extra_compile_args)
extra_link_args = ['-L/home/glaunay/python3.6/lib']
module1 = Extension('ccmap',
                    libraries = ['m'],
                    include_dirs = ['./include', '/home/glaunay/python3.6/include/python3.6m'],
                    sources = ['ccmapmodule.c', './src/mesh.c', './src/transform_mesh.c', './src/decoygen.c','./src/encode.c'],
                    #extra_compile_args=['-D', 'AS_PYTHON_EXTENSION', '-std=c99', '-pedantic', '-L/home/glaunay/python3.6/lib'])#, '-lpython3.6'])
>>>>>>> 2ea47c47a1a2f14fda4ada7312da52c32ad734ee
                    extra_compile_args=extra_compile_args,
		    extra_link_args=extra_link_args
			#extra_compile_args=['-D', 'DEBUG', '-D', 'AS_PYTHON_EXTENSION', '-std=c99', '-pedantic'])
		)

setup (name = 'ccmapModule',
       version = '1.0',
       description = 'This is the C implementation of the mesh based contact map',
       ext_modules = [module1])
