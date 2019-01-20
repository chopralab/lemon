.. _getting_started:

Getting Started
===============

C++
---

The 'fundamental unit' of a **Lemon** workflow in **C++** is the *lambda*
function. All **Lemon** workflows must accept the arguments of
`chemfiles::Frame` and `const std::string&`. An example is given below:

.. code-block:: c++

    auto worker = [](chemfiles::Frame complex, const std::string& pdbid) {
        // do something cool with the complex pdbid and store to my_result
        return my_result;
    };

Note that my_result can be of any type! There is a section on custom, user
defined types later on in the documentation. Once created, the *lambda*
function is passed to the `launch` function.

.. code-block:: c++

    lemon::Options o(argc, argv);
    lemon::launch<lemon::print_combine>(o, worker, std::cout);

The meaning of `print_combine` is to print the results of the workflow using
the `<< operator` and `std::cout`. This may need to change if the user wishes
for a different result or `return`s a custom type from their workflow.

Python
------

The 'fundamental unit' of a **Lemon** workflow in **Python** is the `Workflow`
class provided in the `lemon` module. This class has two functions that should
be `override`n to provide optimal benefit to the author of a workflow. The
first is the `worker` function, taking two additional arguments (other than
`self`). The first is a `chemfiles::Frame` and the second is a string with the
**PDB** id code. Note that both **Python 2** and **Python 3** are supported,
but one must be careful to select the proper version when **Lemon** is build.

.. code-block:: python

    from lemon import *
    class MyWorkflow(Workflow):
        def worker(self, frame,pdbid):
            return get(frame.topology().residue(1), "chainname").get().as_string() + '\n'
        def finalize(self):
            print('Done!\n')

Unlike the **C++** version of the code, the `worker` method *must* return a
string or have no return value. This method is called for all structures in the
**PDB**. The `finalize` method is called when the workflow completes. A
**Python** workflow is invoked by:

.. code-block:: bash

    lemon_python -p myscript.py -w /path/to/mmtf/sequence/files
