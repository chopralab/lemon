.. _count_residues:

Example Programs for counting residues
======================================

A simple start
--------------

In our first example, we will investigate a program that counts the number of
times all residues accur within the PDB.  In this example, we will handle
threading ourselves by creating a map object which holds the residue keys using
the current thread as the key.

.. literalinclude:: ../../progs/count/residues.cpp
   :language: cpp
   :lines: 6-21
   :dedent: 4

Looking at a residue property
-----------------------------

We can also count the number of bioassemblies in a complex:

.. literalinclude:: ../../progs/count/bioassemblies.cpp
   :language: cpp
   :lines: 9-19
   :dedent: 4

or the number of alternative locations in a given complex:

.. literalinclude:: ../../progs/count/altloc.cpp
   :language: cpp
   :lines: 9-19
   :dedent: 4
