.. _prune:

Prune
=====

After residues have been selected, one may wish to remove some residues if they
do not fit a given criterion. The functions below allow one allow to do so.

All functions are availible in both **C++** and **Python**. In **Python**, the
function names should be prefixed with `prune_` instead of the namespace
resolution.

Provided functions
------------------

.. doxygennamespace:: lemon::prune

Example
-------

C++
~~~

The following example demonstrates how to remove cofactors and other 'common'
residues from a selection.

.. literalinclude:: ../../progs/count/small_molecules.cpp
   :language: cpp
   :lines: 9-26
   :dedent: 4

This example extends the previous one to show how one can find only the small-
molecules interacting with a **Heme** group.

.. literalinclude:: ../../progs/interactions/hem_small_molecules.cpp
    :language: cpp
    :lines: 13-31
    :dedent: 4

Python
~~~~~~

These examples are availible in python as:

.. literalinclude:: ../../lang/python/tests/small_molecule.py
    :language: python
    :lines: 1-14

.. literalinclude:: ../../lang/python/tests/small_heme.py
    :language: python
    :lines: 1-26
