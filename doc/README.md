# How to make the documentation:

### HTMD documentation:
The API documentation is generated automatically from the docstrings in the .py
source code files using Sphinx. Any modifications to the docstings will be made
public once the documentation is built and published.

### Creating API sections:
Sphinx automatically generates documentation but it looks too computer generated.
To make it more humanly readable and also to create summary pages, we use reStructuredText. 
Any *.rst file in the htmd/doc/source which does not start with "htmd." is manually written 
and can be modified. RST files starting with "htmd." are produced by Sphinx and will be over-
written on the next build. No need to read or modify those files.

The main documentation tree starts in **htmd/doc/source/index.rst**
Lines under .. toctree:: which are pure text, refer to a RST file and use the 
title defined in that file (1st file line) as the link text when produced in HTML.
Since we cannot modify file titles in Sphinx produced RST files, in some cases
we provide a link text before the link file to make it prettier.

For example 
```
.. toctree::
   :maxdepth: 2

   tutorials  <--------- This is a link to the tutorials.rst file. It will use the file title as link text
   HTMD Protocols <htmd.protocols.rst>   <-------- This is a link to parametrisation.rst file with link text "HTMD Protocols"
```

### Building and publishing:
Use the Makefile in htmd/doc/ to build and publish the documentation with the 
following commands:
```
make clean
make html
```
This produces the htmd in the htmd/doc/build/ folder which can be inspected before
pushing it to the public server.
```
make publish
```
then pushes the new documentation to the www.htmd.org server.

### Static pages
Some pages are written in static html. These can be found in doc/static/ and need to be
modified manually.

### Modifying tutorials:
Tutorials are located in htmd/tutorials/ in the form of Jupyter notebooks.
Change directory to the tutorials folder and start a jupyter notebook server 
from the command line by calling `jupyter notebook`
Now open the notebook and modify anything necessary and then save it.

The Makefile will take any *.ipynb file in htmd/tutorials/ and convert it to an RST
file in htmd/doc/source/tutorials/ folder. This folder is automatically generated
and off-limits for human modification.

To link your tutorial to the htmd website, you need to add it to the htmd/doc/source/tutorials.rst
file and then build.
