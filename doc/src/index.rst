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

C++ header library and example programs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Lemon**'s source code is availible under the BSD license and located on
GitHub_. **Lemon** can be obtained using the following commands in a UNIX-like
environment. To build the **C++** side of **Lemon**, you need a C++11 compiler,
the CMake_ build system.  Note ASYNC features require a C++14 compiler.

For Python support, please install the *Python interpreter* and C development
libraries for the version of Python you wish to use and sure that this version
of **Python** is the default version used on the command-line.

.. _GitHub: http://github.com/chopralab/lemon
.. _CMake: https://cmake.org

.. code-block:: bash

    git clone https://github.com/chopralab/lemon.git
    cd lemon
    mkdir build
    cd build
    # The LEMON_BUILD_PYTHON and LEMON_TEST_ASYNC variables are optional and
    # OFF by default
    cmake .. -DCMAKE_BUILD_TYPE=Release # Give additional build arguments here
    cmake --build . -- #(additional build arguments like -j2 or /m:2)

Prebuilt **Python** module for `pip`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A prebuilt version of **Lemon** is availible on the *Python Package Index*
(PyPI_). It can be installed for recent versions of **Python** (v3.5+). For
Linux platforms, v2.7 is supported as well for legacy reasons. Use the following
command to install it.

.. code-block:: bash

    python3 -m pip install candiy-lemon

.. _PyPI: https://pypi.org/project/candiy-lemon/

Note that this package uses the synchronous threading model and is compiled
using a relatively old version of GCC to meet the requirements of the PyPI
service. Therefore, some users may wish to install the package themselves.
See the next section for details.

Building the **Python** module manually
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Users interested in better compiler options can opt to build the **Lemon**
module themselves. This is the only option for supporting older versions of
**Python** and asynchronous threading. Make sure the `python` command
corresponds to the version of **Python** you wish to build the module for.

.. code-block:: bash

    git clone https://github.com/chopralab/lemon.git
    cd lemon
    python setup.py bdist_wheel --build-type MinSizeRel -- \
      -DPYTHON_EXECUTABLE:FILEPATH=`which python` \
      # Other cmake options here
    python -m pip install *.whl # The name of the wheel file is dependant on the python version

Please be sure to check the output of the build process to ensure the correct
compiler, **Python** interpreter, and other build options are selected.

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
    return

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
