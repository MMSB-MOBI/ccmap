# GL basic makefile
# Add debug and valgrind compiler option

GCC_FLAGS = -std=c99 -Wall -I./include -lm -pedantic 
DEBUG=-g -D DEBUG -O0
VALGRIND_ONLY=-g -O0;
all: ccmap string path

debug: GCC_FLAGS += $(DEBUG)
debug : path

valgrind: GCC_FLAGS += $(VALGRIND_ONLY)
valgrind: path

#src/mesh_default.c deleted
ccmap: src/ccmap.c
	gcc -o bin/$@ $<  src/mesh_default.c src/cell_crawler.c  src/molecular_object.c src/pdb_coordinates.c  src/decoygen.c src/encode.c src/mesh.c src/parameters.c src/transform_mesh.c src/fibonacci.c src/miscellaneous.c src/sasa.c src/atom_mapper.c $(GCC_FLAGS)

string: src/string_test.c
	gcc -o bin/$@ $< $(GCC_FLAGS) 

path: src/path_main.c
	gcc -o bin/$@ $<  src/mesh_default.c src/molecular_object.c src/pdb_coordinates.c\
					  src/my_string.c src/fibonacci.c src/atom_mapper.c src/miscellaneous.c\
					  src/mesh.c src/decoygen.c src/cell_crawler.c src/encode.c src/sasa.c\
					  src/pathfinder.c src/convex.c $(GCC_FLAGS)


clean: 
	rm ./bin/*
