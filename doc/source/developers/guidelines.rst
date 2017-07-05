##############################
Guidelines for HTMD developers
##############################

Coding Style Tips
=================

Note: some style is at variance with Python PIP recommendations.
 
* Class names start with capital letters
* Method names start with lower-case, then ``camelCase``
* Same for generic function, e.g. ``testMe()``
* Modules should be nouns
* Methods and functions should be verbs
* Make the main a test case, where possible.
* Use namespaces instead of composite name, e.g. ``charmm.build()`` instead of ``charmmBuild()`` when possible
* Try to keep single names when possible, so there is no need for camelCase


Using modules inside HTMD
=========================

Creating efficient code is different from generating scripts for personal use. As such, when developing code for HTMD,
be as minimal as possible when importing modules. As examples, please do _not_ use any of these types of imports inside
HTMD code:

.. code:: python

    import htmd
    from htmd import *

Furthermore, from ``from <module> import *`` should _never_ be used, as it pollutes the namespace and can shadow same-name
functionalities. Try as much as possible to only import the function/class you specifically need instead of importing an
entire module, unless one wants to use that module heavily on that implementation. Keep a simple
``htmd.modulename.submodulename`` structure. So file names if they are not meant to be modules (e.g. util.py) should be
imported in the upper module namespace.

Do not pollute the module and submodule names. Changes to the modules structure requires consensus and approval, as well
as documentation creation.


Using Docstrings and Doctests
=============================

See this template:

.. code:: python

    def home(dataDir=None, libDir=False):
        """Return the pathname of the HTMD root directory (or a data subdirectory).

        Parameters
        ----------
        dataDir : str
            If not None, return the path to a specific data directory
        libDir : bool
            If True, return path to the lib directory

        Returns
        -------
        dir : str
            The directory

        Example
        -------
        >>> htmd.home()                                 # doctest: +ELLIPSIS
        '.../htmd'
        >>> htmd.home(dataDir="dhfr")                   # doctest: +ELLIPSIS
        '.../data/dhfr'
        >>> os.path.join(htmd.home(dataDir="dhfr"),"dhfr.pdb")  # doctest: +ELLIPSIS
        '.../data/dhfr/dhfr.pdb'
        """

Docstrings can be test cases (as above). This is convenient because you have four things in one place:

#. the test case
#. the expected result
#. an example
#. the rest of the documentation

It's sufficient to add this in the main:

.. code:: python

    if __name__ == "__main__":
        import doctest

        failure_count, _ = doctest.testmod()
        if failure_count != 0:
            raise Exception('Doctests failed')


The ``doctest: +ELLIPSIS`` comment on the docstring indicates that match with ``...`` is flexible.
Other possibly useful directives are ``SKIP`` and ``NORMALIZE_WHITESPACE``.

One can also:

- run tests placed in external files with ``doctest.testfile('doctest_in_help.rst')``
- test a different module with ``doctest.testmod(doctest_simple)``
