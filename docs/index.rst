.. Example documentation master file, created by
   sphinx-quickstart on Sat Sep 23 20:35:12 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

The Jessie Duncan HSR Project documentation! |:fox:|
====================================================

Welcome! This is the online documentation for the `jdhsrpy <https://github.com/KriSun95/jdhsrpy>`__ code repository.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   setting_up_nustar_data
   other_processing
   other_software
   publications
   auto_examples/index
   modules

Installation
------------

To install the code base, I've decided to us the ``uv`` package manager that can be found `here <https://docs.astral.sh/uv/>`__.

First, clone the repository.

.. code-block:: bash
   :caption: cloning

   git clone https://github.com/KriSun95/jdhsrpy.git

Then, in the top-level of the downloaded repository, the one with ``jdhsrpy`` folder and the ``pyproject.toml`` file, run the following:

.. code-block:: bash
   :caption: installing

   uv pip install .

Once the installation is complete, the code-base can be imported and used Python code.

.. code-block:: python3
   :caption: jdhsrpy

   import jdhsrpy

The final step is to have fun!

Building the documentation
--------------------------

The documentation for the code-base is located on this website, but let's describe how a user might get a local copy.

After installing the repository, you will be able to build build the website documentation easily.

The default lines the user need to run in the same top-level directory descitbed above is to first generate the API documentation for the code-base. In the terminal, run:

.. code-block:: bash
   :caption: api-docs

   sphinx-apidoc -o docs jdhsrpy/

This will create a ``modules.rst`` file in the ``docs`` folder that uses the code doc-strings to build the code descriptions.

Next, to run all the examples and build the erest of the website, in the terminal, run:

.. code-block:: bash
   :caption: local-website-docs

   sphinx-build docs _build

This runs the gallery examples and creates the ``docs/sg_execution_times.rst`` file and ``docs/auto_examples`` folder. It also creates the ``_build`` directory.

The ``_build`` directory is where the full documentation lives. In this folder, the user can open any ``html`` file and then navigate the local website code but it is convention to open ``index.html``.

Rebuilding the documentation
****************************

If the documentation has been built locally but the user wants to build it again, it is annoying but there might be some manual steps needing done.

To save time, not all the code is completely re-run so we should clear all the previous documentation first before re-building.

Therefore, it is safest to just delete the ``docs/sg_execution_times.rst`` file and the ``docs/auto_examples`` and ``_build`` folders before running the ``sphinx`` lines again.

Saving to a custom location
***************************

If you would rather have the documentation website somewhere else, we only need to change the final argument which is currently ``_build``.

For the ``sphinx-apidoc`` command, the two inputs after ``-o`` will have to be these folders. A user can, of course, run this command from anywhere on their machine as long as the inputs still point to those folders, they would just need to include the whole path.

For the ``sphinx-build`` command, the first input still needs to point to the ``docs`` folder wherever that is, but the second input can be anywhere on the user's machine.

For example, I could savee the documentation in website form in my ``Downloads`` folder, in a folder called ``website-docs``, with the following:

.. code-block:: bash
   :caption: docs-location

   sphinx-build docs ~/Downloads/website-docs

This assumes I'm running the command in that same top-level directory to be referencing the first input as ``docs`` but if I was somewhere else on my machine I would just change that to ``<path-to>/docs`` where ``<path-to>`` is the directory location of ``docs``.
