**Lemon**
=========

*a modern C++ library for mining structural information from the PDB*

Introduction
------------

The MacroMolecular Transmission Format (MMTF) is a highly efficient file Format
for the storage and processing of biomolecular data. As of January 2019, the
entire Protein Data Bank (PDB) can be stored in under 10 Gb when using this
format and compressed using the Gzip algorithm. As a result, a copy of the PDB
is availible for here_.

While packages are availible for decoding MMTF structures, these packages do
not provide suitable working environments for exploring and manipulating the
contents of an MMTF file. Further, few packages exist to run PDB-wide
calculations easily, efficiently, and using parallel environments provided by
traditional high-throughput computing clusters. Lemon solves this issue by
providing both a C++-based API for reading PDB Hadoop sequence files and a
parallel environment for launching PDB-wide calculations.
We have deemed these PDB-wide calculations workflows and the remainder of this
document details how to create and run these workflows in both C++ and Python
environments.

.. _here: https://mmtf.rcsb.org/download.html

Obtaining Lemon
---------------

**Lemon**'s source code is availible under the BSD license and located at
github_. **Lemon** can be obtained using the following commands in a UNIX-like
environment. To build the **C++** side of **Lemon**, you need a C++11 compiler,
the CMake_ build system, and the `Boost C++ Libraries`_.  Note ASYNC features
require a C++14 compiler.  For Python support, please install the *Python
interpreter* and C development libraries for the version of Python you wish to
use and sure that this version of Python is the default version used on the
command-line.

.. _github: http://github.com/chopralab/lemon
.. _CMake: https://cmake.org
.. _`Boost C++ Libraries`: https://www.boost.org

.. code-block:: bash

    git clone https://github.com/chopralab/lemon.git
    cd lemon
    mkdir build
    cd build
    # The LEMON_BUILD_PYTHON and LEMON_TEST_ASYNC variables are optional and
    # OFF by default
    cmake .. -DCMAKE_BUILD_TYPE=Release -DLEMON_BUILD_PYTHON=ON -DLEMON_TEST_ASYNC=ON
    make -j 2

Developing a **Lemon** workflow
-------------------------------

Lemon supports workflow development in both C++ or Python. We have documented
for both languages is presented side-by-side as to avoid repetition of the same
text and examples.

The underlying philosophy of the Lemon API is that of functional programming.
Therefore, many functions operate on existing structures that the user
initializes and passes to the operations of interest.

.. toctree::
    :maxdepth: 2

    getting_started
    constants
    select
    prune
    count
    separate

Miscellaneous Functionality
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
    :maxdepth: 2

    geometry
    xscore
    tmalign

Invoking a **Lemon** program
----------------------------

.. toctree::
    :maxdepth: 2

    options
    parallel

Example workflows
-----------------

.. toctree::
    :maxdepth: 1

    count_residues
    small_molecules
