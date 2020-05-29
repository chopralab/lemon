.. tmalign:

TMalign
=======

Some users may wish to write **Lemon** workflows that involve the alignment of
a potential protein to a reference protein. Potential applications include the
creation of a benchmarking set of all proteins with a given fold or motif. To
address this, we've included the *TMalign* algorithm within **Lemon**.

**Note:** You will need to include the additional file: `<lemon/tmalign.hpp>`
to have access to these features.

The TMscore function for alignment
----------------------------------

.. doxygenstruct:: lemon::tmalign::TMResult

.. doxygenfunction:: lemon::tmalign::TMscore

Example
-------

C++
~~~

.. literalinclude:: ../../progs/misc/tmscore_all.cpp
   :language: cpp
   :lines: 8-29
   :dedent: 4

Python
~~~~~~

.. literalinclude:: ../../lang/python/tests/tmalign.py
   :language: python
   :lines: 1-16
