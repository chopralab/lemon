## Lemon:  Tool for developing benchmarking sets from the entire PDB in minutes 

[![Build Status](https://travis-ci.org/chopralab/lemon.svg?branch=master)](https://travis-ci.org/chopralab/lemon)

[![Coverage Status](https://coveralls.io/repos/github/chopralab/lemon/badge.svg?branch=master)](https://coveralls.io/github/chopralab/lemon?branch=master)

### What is Lemon's purpose?

**Lemon** is a tool for developing benchmarking sets for structural biology software dealing with 3D Macromolecules.  It is designed to be fast and flexible, allowing users to quickly query the 3D features of a given collection of 3D structures.  To do so, the user writes a C++11 Lambda function which is applied to all the structures selected by the user.

Due to the incredibly fast parsing speed of the MMTF format, Lemon uses this format by default.  This helps **Lemon** query the entire Protein Data Bank under 25 minutes on an 8 core machine.

### How do I obtain Lemon?

**Lemon** is developed to have as few dependencies as possible. You only need a recent C++ compiler which supports C++ and a copy of the Boost Filesystem library. All other dependencies are installed for you by the build system.

```bash
git clone https://gitlab.com/chopralab/lemon.git

cd lemon

mkdir build

cd build

cmake .. -DCMAKE_BUILD_TYPE=Release

make -j 2

```

### How does one use Lemon?

The Protein Data Bank is used to test **Lemon**'s capabilites and is the source of the majority of structural biology benchmarking sets.  Therefore we have included a script to download the entire PDB archive.  It is recommended to use the latest Hadoop sequence files located [here](https://mmtf.rcsb.org/v1.0/hadoopfiles/full.tar).

Currently, the archive takes ~9Gb of space.

**Lemon** can use create your own `idx` files using search queries on [RCSB](https://rcsb.org) and downloading the search result.

To run **Lemon**, select a program. For example, if one wants to query all the small molecules which interact with `SAM`, use the following command:

```bash
tar xf full.tar /dev/shm/
/path/to/lemon/build/progs/count_sam_small_molecules -w /dev/shm/full -n <number of cores>
```

The results are printed to `stdout`.

**Lemon** is &copy; 2018 Chopra Lab and Purdue University. Developed by Jonathan Fine and is  available as open source under the terms of the [BSD License](http://opensource.org/licenses/BSD). 
