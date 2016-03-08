
Python primer for HTMD
======================

by G De Fabritiis

Basics
------

1. Downloaded and install HTMD from www.htmd.org
2. We use python3 syntax.
3. This is a notebook (walk around)

Python imports
--------------

In Python modules have to be imported in order to be used

.. code:: python

    #import basic modules from interactive use
    from htmd import *    
    #Optional: set default plot to be inline
    %matplotlib inline            


.. parsed-literal::

    HTMD. 26-27 November 2015 HTMD workshop in Barcelona (fully booked)
    
    You are on the latest HTMD version (0.1.11).


.. code:: python

    a = 3+2
    print('a is', a)


.. parsed-literal::

    a is 5


Notebooks
---------

-  python is the interpreter
-  ipython is an interactive interpreter
-  A notebook like this one is a grapical interactive environment.
-  To navigate the filesystem (% indicates a notebook command)

   .. code:: python

       %ls
       %mkdir test

Interacting with the filesystem
-------------------------------

Normal filesystem commands work from the interactive session but in a
program the functional version should be used

.. code:: python

    os.chdir('.') 
    os.listdir('.')
    glob('htmd-*')




.. parsed-literal::

    ['htmd-python-primer.ipynb']



Strings
-------

Strings are identified by single or double apices like in

.. code:: python

    a='newtest' 
    b="newtest"

Relational operators
--------------------

Relational operators are ``==,<,>,<=,>=,!=`` and can be connected with
``and,or,not``

.. code:: python

    a==b




.. parsed-literal::

    True



Lists
-----

Lists of different objects are created using square brackets

.. code:: python

    a = [1, 2, 3, 'ba', 5, 6, 7, 8, 9]

Help
----

The command help provides help for a given function

.. code:: python

    help(print)


.. parsed-literal::

    Help on built-in function print in module builtins:
    
    print(...)
        print(value, ..., sep=' ', end='\n', file=sys.stdout, flush=False)
        
        Prints the values to a stream, or to sys.stdout by default.
        Optional keyword arguments:
        file:  a file-like object (stream); defaults to the current sys.stdout.
        sep:   string inserted between values, default a space.
        end:   string appended after the last value, default a newline.
        flush: whether to forcibly flush the stream.
    


.. code:: python

    #shift tab
    print(


::


      File "<ipython-input-8-30fdca73e0d6>", line 2
        print(
              ^
    SyntaxError: unexpected EOF while parsing



Arrays
------

Arrays are best created using numpy array and plots using matplotlib.

.. code:: python

    m = np.array([[ 11., 12, 13, 14 ],
               [ 21, 22, 23, 24 ],
               [ 31, 32, 33, 34 ]])
    print(m)


where a is 4-by-4 matrix of double numbers, due to the fact than ``11.``
is a real number.

Array multiplication
--------------------

Arrays can be multiplied easily element by element

.. code:: python

    b = 3.0 * m

A scalar multiplication is applied to each element of the array.

Sequences
---------

An array sequence can be created with ``arange``. Array indexing starts
at 0.

.. code:: python

    a=np.arange(0,4)
    print(a)

Slicing
-------

.. code:: python

    a[1:] #from 1 till the end, (starts at zero)
    m[:,0] = 99 # first column
    m[-2:,] #indexing backwards is possible
    a[0:4] # from 0 to 3

Some matrix operations
----------------------

.. code:: python

    np.multiply(m,a)

.. code:: python

    m.transpose()

.. code:: python

    np.inner(m,a)

Other operations
----------------

.. code:: python

    np.concatenate((a,a))

.. code:: python

    m.sum(axis=0)

Array properties
----------------

Basic information about the array can be accessed by

.. code:: python

    print(a.shape)
    print(a.size)
    print(a.ndim)

Simple plots
------------

.. code:: python

    from matplotlib.pyplot import plot,title,xlabel,ylabel,grid,show,savefig
    x = np.arange(0,50,.5)
    y = np.sin(x/3) - np.cos(x/5)
    plot(x,y, '.-k')
    plot(x,np.sin(x/2) - np.cos(x/5),'.r')
    title('A simple double plot')
    xlabel('variable 1'), ylabel('variable 2')
    grid(True)
    show()
    savefig('/tmp/foo.eps')

Conditionals
------------

.. code:: python

    if 1>2:
        a=100
        print('2')

as you probably know in python spaces are important, use indentation to
define a scope after the :

Loops
-----

.. code:: python

    for i in np.arange(1,6):
        print(i)

Functions
---------

Subroutines are defined using ``def``

.. code:: python

    def test(a,b=1,c=3,d=1):
        return a*b*c
    
    test(1,c=5)

Argument passing
----------------

-  Python represents all its data as objects.
-  Variables are just names.
-  Some objects are mutable, some immutable.
-  Immutables are: int, float, complex, str, tuples, bytes, frozensets
-  Mutables are: list, byte array, sets, dict, classes

Identity of an object
---------------------

With ``id`` is possible to check the unique

.. code:: python

    n = 1 # integer immutable
    id(n)

.. code:: python

    n += 1
    id(n) # new object

.. code:: python

    m = [1] #list mutable
    id(m)

.. code:: python

    m.append(2)
    id(m) #same object

Argument passing
----------------

-  Passing an argument to a function is like creating a new name to the
   object.
-  If it is mutable then any change inside the function will affect the
   object outside.
-  If it is immutable and the function changes it, then python creates
   another object inside the function scope, so nothing changes outside
   of the function.

Debugging code in notebooks
---------------------------

Debug an error and see what's wrong by inspecting variables

.. code:: python

    def broken_function(b):
        a = 1
        print(a,b.xxx())
    c = 4
    broken_function(c)

.. code:: python

    %debug  

 HTMD documentation
-------------------

-  Embedded inside the code
-  Available online at www.htmd.org




Exercises
---------

1. Sum the first 50 numbers with a for loop
2. Do the same using numpy arrays
3. Write a function that sets a value of its argument (that is an
   integer)
4. The same but for a numpy array
