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

~~~ 


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

~~~


Doctests
--------

Docstrings can be test cases (as above). The "ELLIPSIS" line indicates that 
match with dots is flexible.   

It's sufficient to add this in the main...
 
~~~
if __name__ == "__main__":
    import doctest
    doctest.testmod()
~~~

