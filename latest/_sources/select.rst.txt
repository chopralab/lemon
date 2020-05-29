.. _select:

Select
======

Most **Lemon** workflows start by selecting one or multiple sets of residues
and performing operations on these residues. The functions below are availible
in both **C++** and **Python**, but there are a few implementation differences
users should note.

First, the **C++** version of all functions are implemented as templates,
allowing for the user to use their prefered container for storing the resulting
residues. Second, each function has two overloads. One takes a container
initialized by the user as an argument and subsequently populates it using the
container's `insert` method. No template arguments are required as these can be
deduced at compile-time. An other overload initializes and populates a new
container specified by a template argument (default is a `std::list<size_t>`).
The choice of correct container is left to the user.

In **Python**, both overloads are availible as well. However, due to
restrictions imposed by **Python** generics, the user must use the `ResidueIDs`
container. All **Python** functions have been prepended with *select_*.


Provided selectors
------------------

.. doxygennamespace:: lemon::select

Example
-------

C++
~~~

.. literalinclude:: ../../progs/count/metal_ions.cpp
   :language: cpp
   :lines: 9-19
   :dedent: 4

Python
~~~~~~

.. literalinclude:: ../../lang/python/tests/metal_ions.py
    :language: python
    :lines: 1-14
