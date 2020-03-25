# GL basic makefile
# Add debug and valgrind compiler option

GCC_FLAGS = -std=c99 -Wall -I./include -lm#-pedantic 
DEBUG=-g -D DEBUG -O0

all: ccmap string

debug: GCC_FLAGS += $(DEBUG)
debug : ccmap

ccmap: src/ccmap.c
	gcc -o bin/$@ $<  src/cell_crawler.c src/default.c src/molecular_object.c \
								  src/pdb_coordinates.c  src/decoygen.c src/encode.c src/mesh.c \
								  src/parameters.c       src/transform_mesh.c $(GCC_FLAGS)

string: src/string_test.c
	gcc -o bin/$@ $< $(GCC_FLAGS) src/default.c
     
clean: 
	rm ./bin/*
