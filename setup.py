#from distutils.core import setup, Extension
from setuptools import setup, find_packages, Extension
import os
import sysconfig

extra_compile_args = sysconfig.get_config_var('CFLAGS').split()
extra_compile_args += ['-D', 'AS_PYTHON_EXTENSION', '-std=c99', '-pedantic'] #'-D', 'DEBUG', '-D', 'PYMEM_CHECK']
core = Extension('ccmap',
              libraries = ['m'],
              include_dirs = ['ccmap/include'],
              sources = ['ccmap/ccmapmodule.c', 'ccmap/src/ccmapmodule_utils.c', 'ccmap/src/ccmapmodule_allocation.c','ccmap/src/decoygen.c', 'ccmap/src/encode.c', 'ccmap/src/molecular_object.c', 'ccmap/src/pdb_coordinates.c', 'ccmap/src/cell_crawler.c', 'ccmap/src/mesh.c', 'ccmap/src/transform_mesh.c', 'ccmap/src/mesh_default.c' ],
              extra_compile_args=extra_compile_args
		)
setup (name = 'ccmap',
       version = '2.1.2',
       author = 'G.Launay', 
	author_email='pitooon@gmail.com',
       description = 'A C implementation of a mesh based atomic pairwise distance computating engine, with docking pose generation capabilities',
       ext_modules = [core],
       long_description_content_type='text/markdown',
       long_description=open('README.md').read(),
       url="https://github.com/MMSB-MOBI/ccmap",
       keywords="protein docking bioinformatics structure",
       classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Development Status :: 5 - Production/Stable" 
       ],
       )


