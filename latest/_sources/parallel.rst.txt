.. _invoke-workflow:

Launch a **Lemon** workflow
===========================

C++ Workflows
-------------

Since **Lemon** is a *header-only* library, users of the C++ API have complete
control over how they wish to launch their workflows. The `lemon::Options`
class allows users to read command line options in a consistant manner. See the
previous section for a more indepth discussion of this class. The following
function is used to launch a lemon workflow.

.. doxygenfunction:: lemon::launch

Submitting **Lemon** jobs
~~~~~~~~~~~~~~~~~~~~~~~~~

If a **Lemon** workflow does not require any custom features (IE custom options
passed to the `main` function from the shell), users can use the script called
`launch_lemon.pbs` to launch their workflow on PBS-based submission systems. An
example of how to use this given below for a sample workflow.

.. code-block:: bash

    wget -N https://mmtf.rcsb.org/v1.0/hadoopfiles/full.tar
    qsub launch_lemon.pbs -v LEMON_PROG=protein_angle

**Note:** The provided script will *not* work on all systems and may need to be
edited to work in a given computing environment.

Python Workflows using lemon_python
-----------------------------------

The program `lemon_python` is secretly a **Lemon** workflow! This workflow
searches for a class called `MyWorkflow` (which must subclass `lemon.Workflow`.
It then instantiates this class and runs the member function `Workflow.worker`
for every entry in the PDB. The member function `Workflow.finalize` is called
when the workflow completes. An example of using this workflow is given below:

.. code-block:: bash

    wget -N https://mmtf.rcsb.org/v1.0/hadoopfiles/full.tar
    tar xf full.tar
    lemon_python -p tmscore.py -w full/

**Note:** The entire python script passed to `lemon_python` is evaluated before
the `MyWorkflow` class is instantiated. Control is not passed back to the
script after the `Worker.finalize` method is called.

Using the Python Interpreter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For users of the `candiy_lemon` PyPI package, one must submit the derived
**Lemon** workflow manually using the `lemon.launch` function. An example is
given below. The arguments to this function are the `lemon.Workflow` daughter
class, the path to the RCSB Hadoop files and the number of cores to use.

.. code-block:: python

    from candiy_lemon.lemon import *
    class MyWorkflow(Workflow):
        def worker(self, entry, pdbid):
            return entry.topology().residue(1).get("chainname").get().as_string() + '\n'
        def finalize(self):
            pass

    work = MyWorkflow()

    launch(work, "full", 2)

Prefiltering the PDB with searches originating on RCSB
------------------------------------------------------

The advanced search features on the RCSB website can be used to prefliter the
PDB. First, one performs a search on the website. Then, they must obtain the
query details by clicking the blue button on the bottom left corner of the
search results. This will result in an XML version of their search being
generated. This must be given as input to the script
`obtain_entries_from_search.pl` provided with **Lemon** which will produce
an entry file that can be passed to a workflow. An example *XML* result and
entry file example are given below:

.. code-block:: xml

    <orgPdbQuery>
        <version>head</version>
        <queryType>org.pdb.query.simple.AdvancedKeywordQuery</queryType>
        <description>Text Search for: hiv</description>
        <queryId>AC7A2BC3</queryId>
        <resultCount>2875</resultCount>
        <runtimeStart>2019-01-28T00:51:01Z</runtimeStart>
        <runtimeMilliseconds>1686</runtimeMilliseconds>
        <keywords>HIV</keywords>
    </orgPdbQuery>

Then

.. code-block:: bash

    perl obtain_entries_from_search.pl hiv_search.xml > hiv_prots.lst
    tar xf full.tar
    ./small_molecules -w full -e hiv_prots.lst

Danger Zone: Internal documentation!
------------------------------------

These functions are internal to **Lemon** and not meant to be used as part of
the external API. They are documented here for the interested reader and future
**Lemon** developers.

.. doxygenfunction:: lemon::run_parallel

.. doxygenclass:: lemon::Hadoop
    :members:

.. doxygenfunction:: lemon::read_hadoop_dir
