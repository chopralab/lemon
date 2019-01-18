.. tmalign:

TMalign
=======

Some users may wish to write **Lemon** workflows that involve the alignment of
a potential protein to a reference protein. Potential applications include the
creation of a benchmarking set of all proteins with a given fold or motif. To
address this, we've included the *TMalign* algorithm within **Lemon**.

The TMscore function for alignment
----------------------------------

.. doxygenstruct:: lemon::tmalign::TMResult

.. doxygenfunction:: lemon::tmalign::TMscore

Example
-------

.. literalinclude:: ../../progs/misc/tmscore_all.cpp
   :language: cpp
   :lines: 7-26
   :dedent: 4
