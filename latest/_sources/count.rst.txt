.. _calculation-count:

Counting operations
===================

These operations are provided so that the user can easily count the number of
times a given property occurs. All operations are availible in **C++** and
**Python**. **Python** functions are prefixed with `count_` instead of using
the namespace resolution.

Provided selectors
------------------

.. doxygennamespace:: lemon::count

Examples
--------

C++
~~~

.. literalinclude:: ../../progs/count/residues.cpp
   :language: cpp
   :lines: 6-23
   :dedent: 4

Python
~~~~~~

.. literalinclude:: ../../lang/python/tests/residues.py
   :language: python
   :lines: 1-12
