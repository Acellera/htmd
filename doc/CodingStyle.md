Guidelines for developers
=========================

Note: some style is at variance with Python PIP recommendations.
 
* Class names start with capital letters
* Method names start with lower-case, then camelCase
* Same for generic function,e.g.   testMe()
* Make the main a test case, where possible. 
* Name methods as verbs.
* Use namespaces instead of composite name, e.g. charmm.build() instead of charmmBuild() when possible



Docstrings
----------

See this template:


```python
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

```


Doctests
--------

Docstrings can be test cases (as above). This is convenient because you have four 
things in one place: 

 1. the test case
 2. the expected result
 3. an example
 4. the rest of the documentation

It's sufficient to add this in the main...
 
```python
if __name__ == "__main__":
    import doctest
    # ... add any other global / import you want tests to see
    doctest.testmod()
```

The "ELLIPSIS" line indicates that match with dots is flexible.   
Other possibly useful directives are SKIP and NORMALIZE_WHITESPACE. You can also

  * run tests placed in external files with `doctest.testfile('doctest_in_help.rst')`
  * test a different module with `doctest.testmod(doctest_simple)`

