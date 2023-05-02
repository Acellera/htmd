Installation
============

If you are a non-profit entity, you can download a free version of HTMD which does not have all functionalities of the Pro-version.

The Pro version requires a license, which is linked to acceptance of the end user License agreement (EULA).

If you do not have a license, please contact us. If you already have a license, please follow the instructions here below.

Installing HTMD
---------------

If you do not have conda installed please first install Miniconda by downloading and running the latest installer from the following URL https://docs.conda.io/en/latest/miniconda.html

After installing conda create a new environment for HTMD and activate it.

.. code-block:: bash
    :linenos:

    conda create -n htmd
    conda activate htmd

With the ``htmd`` environment active (you should see ``(htmd)`` on the left) run the following command to install the HTMD python package

.. code-block:: bash
    :linenos:

    conda install mamba python=3.10 -c conda-forge
    mamba install htmd -c acellera -c conda-forge

Before using HTMD, please read the EULA and register using the command

.. code-block:: bash
    :linenos:

    htmd_register

If you don't accept the EULA and register at this point, HTMD will prompt for registration when the Python module is imported 