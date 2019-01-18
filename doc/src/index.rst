**Lemon**
=========

*a modern C++ library for mining structural information from the PDB*

Introduction
------------

The MacroMolecular Transmission Format (MMTF) is a highly efficient file Format
for the storage and processing of biomolecular data.  As January 2019,
the entire Protein Data Bank (**PDB**) can be stored in under *10 Gb* when
using this format and compressed using the **Gzip** algorithm.  As a result, a
copy of the **PDB** is availible for :download:`here <https://mmtf.rcsb.org/download.html>`.

While packages are availible for decoding MMTF structures, these packages do
not provide suitable working environments for exploring and manipulating the
contents of an MMTF file. Further, few packages exist to run *PDB-wide*
calculates easily, efficiently, and using parallel environments provided by
traditional high-throughput computing clusters. **Lemon** solves this issue by
providing both a C++-based API for reading **PDB** Hadoop sequence files, and
a parallel environment for launching PDB-wide calculations.

We have deemed these PDB-wide calculations *workflows* and the remainder of
this documentation details how to create and run these *workflows* in both
**C++** and **Python** environments.

Developing a **Lemon** workflow
-------------------------------

**Lemon** workflows can be developed in either **C++** or **Python**. The
documentation for both languages is presented *side-by-side* as to avoid
repatition of the same text and examples.

The underlying philosophy of the **Lemon** API is that of functional
programming. Therefore, many functions operate on existing structures that the
user initializes and passes to the operations they are interested in.

.. toctree::
    :maxdepth: 2

    getting_started
    constants
    select
    prune
    count
    separate

Miscellaneous Functionality
---------------------------

.. toctree::
    :maxdepth: 2

    xscore
    tmalign

Invoking a **Lemon** program
----------------------------

.. toctree::
    :maxdepth: 2

    options
    hadoop
    parallel

Example workflows
-----------------

.. toctree::
    :maxdepth: 1

    count_residues
    small_molecules
