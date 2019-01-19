.. xscore:

XScore/Vina
===========

When developing a benchmarking set, a user may wish to select complexes with a
certain set of interactions.  We provide the scoring function used by the
AutoDOCK vina software package to help users select desirable interactions.

.. doxygenstruct:: lemon::xscore::VinaScore

.. doxygenfunction:: lemon::xscore::vina_score

Example
-------

C++
~~~

.. literalinclude:: ../../progs/misc/vinascore_all.cpp
   :language: cpp
   :lines: 9-49
   :dedent: 4

Python
~~~~~~

.. literalinclude:: ../../lang/tests/vina.py
    :language: python
    :lines: 1-42
