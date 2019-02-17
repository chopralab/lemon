.. _getting_started:

Getting Started
===============

Currently, **Lemon** supports development of *workflows* in both **C++** and
**Python**. Note that **Python** development is split into two modes which are
described in the documentation below.

C++
---

Creating a new *workflow*
~~~~~~~~~~~~~~~~~~~~~~~~~

**Lemon** is a header library that requires the use of the
*chemfiles* library (which is compiled). This means that **Lemon** can be used
in *chemfiles* projects by simply copying the include directory. For simple
*workflows*, one can add them to the **Lemon** build tree by creating a *cpp*
file in one of the *progs* subdirectories (you will need to edit the
*CMakeLists.txt* file if you create a new subdirectory). When CMake is run, it
should automatically find your *cpp* file and add it to the build tree.

Developing the *workflow*
~~~~~~~~~~~~~~~~~~~~~~~~~

The 'fundamental unit' of a **Lemon** workflow in **C++** is the *lambda*
function. All **Lemon** workflows must accept the arguments of
`chemfiles::Frame` and `const std::string&`. An example is given below:

.. code-block:: c++

    auto worker = [](chemfiles::Frame entry, const std::string& pdbid) {
        // do something cool with the entry, pdbid and store to my_result
        return my_result;
    };

Note that `my_result` can be of any type! There is a section on custom, user
defined types later on in the documentation. Once created, the *lambda*
function is passed to the `launch` function.

.. code-block:: c++

    lemon::Options o(argc, argv);
    auto collector = lemon::print_combine(std::cout);
    lemon::launch(o, worker, collector);

The meaning of `print_combine` is to print the results of the workflow using
the `operator<<` and `std::cout`. This may need to change if the user wishes
for a different result or returns a custom type from their workflow.

Python
------

Development in **Python** can be accomplished in two ways. First, users can use
the provided *pip* package. This is the recommend way to run **Lemon** in an
environment where **Python** development is prefered. Second, users can the
`lemon_python` program which is created when the **C++** code is compiled. This
program is recommend for **C++** developers who like to prototype code in
**Python** before coverting it to a **C++** only version.

Note that there are minor differences between the two versions that are
actively being resolved.

Using the *PyPI* `pip` package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The 'fundamental unit' of a **Lemon** workflow in **Python** is the `Workflow`
class provided in the `candiy_lemon.lemon` module. This class has two functions
that should be `overriden` to provide optimal benefit to the author of a
workflow. The first is the `worker` function, taking two additional arguments
(other than `self`). The first is a `lemon.Frame` and the second is a
string with the **PDB** id code.

.. code-block:: python

    from candiy_lemon import lemon # You must import lemon from the CANDIY module

    min_size = 10

    class MyWorkflow(lemon.Workflow): # Lemon.Workflow must be subclassed
        def worker(self, entry, pdbid):
            # No need to reimport any objects or other items.
            # All variables are carried over
            heme_names = lemon.ResidueNameSet()

            smallm = lemon.select_small_molecules(entry, lemon.small_molecule_types, min_size)

            return pdbid + lemon.count_print_residue_names(entry, smallm) + '\n'

        def finalize(self):
            pass

    a = MyWorkflow()

    lemon.launch(a, "rcsb_hadoop", 6)

Using the `lemon_python` command
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The 'fundamental unit' of a **Lemon** workflow in **Python** is the `Workflow`
class provided in the `lemon` module. This class has two functions that should
be `overriden` to provide optimal benefit to the author of a workflow. The
first is the `worker` function, taking two additional arguments (other than
`self`). The first is a `lemon.Frame` and the second is a string with the
**PDB** id code. Note that both **Python 2** and **Python 3** are supported,
but one must be careful to select the proper version when **Lemon** is built.

.. code-block:: python

    import lemon # Note that we are importing lemon directly
    class MyWorkflow(lemon.Workflow): # You must subclass Workflow
        def worker(self, frame,pdbid):
            # Globals are not availible due to CPython not working properly
            import lemon # You must reimport lemon here do to CPython silliness
            return get(frame.topology().residue(1), "chainname").get().as_string() + '\n'
        def finalize(self):
            print('Done!\n')

Unlike the **C++** version of the code, the `worker` method *must* return a
string or have no return value. This method is called for all structures in the
**PDB**. The `finalize` method is called when the workflow completes. A
**Python** workflow is invoked by:

.. code-block:: bash

    lemon_python -p myscript.py -w /path/to/mmtf/sequence/files
