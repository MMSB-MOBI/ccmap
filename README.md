# A Python package and C Library for fast molecular contact map computation

[Current Version 2.1.3](https://pypi.org/project/ccmap/)

This package was designed as a tool to quickly compute thousands of sets of atomic or residue molecular contacts. The contacts can be evaluated inside a single body or across two bodies. The library scales well, with the support of the native python multithreading.
The module also provides docking poses evaluation by the application of triplets of Euler angles and translation vectors to initial unbound conformations.

## Installing and using the python module

### Installation

Should be as simple as `pip intstall ccmap`. Alternatively you can clone this repo and run `python setup.py install` at the root folder.
Current release was successfully installed through pip on the following combinations of interpreter/platforms.

* python3.8/OSX.10.14.6
* python3.8/Ubuntu LTS

### Usage

From there you can load the package and display its help.

```python
import ccmap
help(ccmap)
```

#### Functions

Four functions are available:

* cmap: computes the contacts of one single/two body molecule
* lcmap: computes the contacts of a list of single/two body molecules
* zmap: computes the contacts between a receptor and a ligand molecule after applying transformations to the ligand coordinates
* lzmap: computes many sets of contacts between a receptor and a ligand molecule, one for each applied ligand transformation

#### Parameters

All module functions take molecular object coordinates as dictionaries, where keys are atoms descriptors and values are lists.

* 'x' : list of float x coordinates
* 'y' : list of float x coordinates
* 'y' : list of float x coordinates
* 'seqRes' : list of strings
* 'chainID' : list of one-letter string
* 'resName' : list of strings
* 'name' : list of strings

#### Additional arguments

##### Contact threshold distance

In Angstrom's unit, its default value is 4.5. It can be redefined by the name parameter `d`.

##### encode : Boolean

If True, contacts are returned as integers. Each integer encoding one pair of atoms/residues positions in contact with this simple formula,

```python
def K2IJ(k, sizeBody1, sizeBody2):
    nCol = sizeBody2 if sizeBody2 else sizeBody1
    return int(k/nCol), k%nCol
```

if False, contacts are returned as strings of JSON Objects

##### atomic : Boolean

If True, compute contact at the atomic level. By default, this if False and the contacts are computed at the residue level.

##### apply : Boolean

If True, the past dictionaries of coordinates will be modified according to Euler/translation parameters.
This is useful to generate single docking conformation.
This argument is only available for the **cmap** function.

##### offsetRec and offsetLig

When working with protein docking data, unbound conformations are often centered to the origin of the coordinates system. Specify the translation vectors for each body with the `offsetRec` and `offsetLig` named arguments. Only available for the **zmap** and **lzmap** functions.

### Working with PDB coordinates files

#### Parsing coordinate data

We usually work with molecules in the PDB format. We can use the [pyproteinsExt](https://pypi.org/search/?q=pyproteinsExt) package to handle the boilerplate. 

```python
import pyproteinsExt
parser = PDB.Parser()
pdbREC = parser.load(file="dummy_A.pdb")
pdbDictREC = pdbREC.atomDictorize
pdbDictREC.keys()
#dict_keys(['x', 'y', 'z', 'seqRes', 'chainID', 'resName', 'name']) ```
```

By convention, following examples will use two molecules names REC(eptor) and LIG(and).

```python
pdbLIG = parser.load(file="dummy_B.pdb")
pdbDictLIG = pdbLIG.atomDictorize
pdbDictLIG.keys()
#dict_keys(['x', 'y', 'z', 'seqRes', 'chainID', 'resName', 'name']) ```
```

## Examples

### Computing single body contact map

#### Computing one map

Setting contact distance of 6.0 and recovering residue-residue contact as an integer list.

```python
ccmap.cmap(pdbDictLIG, d=6.0, encode=True)
```

#### Computing many maps

Using default contact distance and recovering atomic contact maps as JSON object string. The first positional argument specifies a list of bodies to process independently. 

```python
import json
json.load( ccmap.lcmap([ pdbDictLIG, pdbDictREC ], atomic=True) )
```

### Computing two-body contact map

#### Straight computation of one map

The second positional argument of **cmap** is optional and defines the second body.

```python
ccmap.cmap(pdbDictLIG, pdbDictLIG, d=6.0, encode=True)
```

#### Straight computation of many maps

The second positional argument of **lcmap** is an optional list of second bodies. The first two arguments must be of the same size, as the *i*-element of the first will be processed with the *i*-element of the second.

```python
ccmap.lcmap([pdbDictREC_1, ..., pdbDictREC_n], [pdbDictLIG_1, pdbDictLIG_n], d=6.0, encode=True)
```

#### Computation of one map after conformational change

Use the **zmap** function with third and fourth positional arguments respectively specifying the :

* Euler angles triplet
* translation vector

```python
ccmap.zmap(pdbDictREC, pdbDictLIG , (e1, e2, e3), (t1, t2, t3) )
```

Transformations are always applied to the coordinates provided as a second argument, e.g. : `pdbDictLIG`.

#### Computation of many maps after conformational changes

Use the **lzmap** function, arguments are similar but for the Euler angles and translation vectors which must be supplied as lists.

```python
ccmap.lzmap(pdbDictREC, pdbDictLIG , [(e1, e2, e3),], [(t1, t2, t3),] )
```

### Generating docking conformations

The conformations obtained by coordinate transformation can be back mapped to PDB files.
Here, offset vectors `[u1, u2, u3]` and `[v1, v2, v3]` respectively center `pdbDictREC` and `pdbDictLIG` and one transformation defined by the `[e1, e2, e3]` Euler's angles and the `[t1, t2, t3]` translation vector is applied to `pdbDictLIG`. The resulting two-body conformation is finally **applied** to the provided `pdbDictREC` and `pdbDictLIG`. These updated coordinates update the original PDB object for later writing to file.

```python
# Perform computation & alter provided dictionaries
ccmap.zmap( pdbDictREC, pdbDictLIG,
\ [e1, e2, e3], [t1, t2, t3],
\ offsetRec=[u1, u2, u3],
\ offsetLig=[v1, v2, v3],
\ apply=True)
# Update PDB containers from previous examples
pdbREC.setCoordinateFromDictorize(pdbDictREC)
pdbLIG.setCoordinateFromDictorize(pdbDictLIG)
# Dump to coordinate files
with open("new_receptor.pdb", "w") as fp:
    fp.write( str(pdbREC) )
with open("new_ligand.pdb", "w") as fp:
    fp.write( str(pdbLIG) )
```

## Multithreading

The C implementation makes it possible for the ccmap functions to release Python Global Interpreter Lock. Hence, "actual" multithreading can be achieved and performances scale decently with the number of workers. For this benchmark, up to 50000 docking poses were generated and processed for three coordinate sets of increasing number of atoms: 1974([1GL1](https://www.rcsb.org/structure/1GL1)) 3424([1F34](https://www.rcsb.org/structure/1F34)) 10677([2VIS](https://www.rcsb.org/structure/2VIS)).

<figure> <img src="notebook/img/LZMAP_benchmark_1.png" alt="benchmark" /> </figure>

A simple example of a multithread implementation can be found in the provided [script](tests/scripts/threadsTest.py). The `tests` folder allows for the reproduction of the above benchmark.

## Installing and using the C Library

C executable can be generated with the provided makefile. The low-level functions are the same, but the following limitations exist:

* One computation per executable call
* No multithreading.

# Finding the optimal molecular path connecting two atoms
Using a thiner mesh size, it is possible to obtain the shortest path connecting two atoms.
The atoms to connect must be solvent accessible and path search will operate over the surface and solvant accessible cells. The solvent excluded volume is computed for each atom as the sum of of its Van Der Waals radius and the radius of water molecule.
## install
```sh
 git clone -b fibonacci git@github.com:MMSB-MOBI/ccmap.git
 cd ccmap/ccmap
 make pathfinder
```
## Usage

`./bin/pathfinder -x 'ILE:A:1:CA' -y 'LYS:A:13:N' -i ../tests/structures/small_peptide_noH_1model.pdb`

## Effect of parameters on search
The path is guaranteed to be optimal but may take some time to run.
TO COMPLETE
### cell size
### water probe
