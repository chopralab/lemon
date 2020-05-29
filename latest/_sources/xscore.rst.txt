.. xscore:

XScore/Vina
===========

When developing a benchmarking set, a user may wish to select entries with a
certain set of interactions.  We provide the scoring function used by the
AutoDOCK vina software package to help users select desirable interactions.

**Note:** You will need to include the additional file: `<lemon/xscore.hpp>` to
have access to these features.

.. doxygenstruct:: lemon::xscore::VinaScore

.. doxygenfunction:: lemon::xscore::vina_score

Example
-------

C++
~~~

.. literalinclude:: ../../progs/misc/vinascore_all.cpp
   :language: cpp
   :lines: 10-51
   :dedent: 4

Python
~~~~~~

.. literalinclude:: ../../lang/python/tests/vina.py
    :language: python
    :lines: 1-42
