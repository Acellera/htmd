# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import nbformat
from nbconvert import PythonExporter
from glob import glob
import os
import pytest


def test_tutorials(testfolder, tutorials=('ligand-binding-analysis', 'protein-folding-analysis')):
    startdir = os.path.dirname(os.path.realpath(__file__))

    pyexp = PythonExporter()
    pyexp.template_file = os.path.join(startdir, 'simplepython.tpl')

    # For each tutorial noteboook
    for tut in glob(os.path.join(startdir, '*.ipynb')):
        name = os.path.splitext(os.path.basename(tut))[0]

        if name not in tutorials:
            continue

        testsubf = os.path.abspath(os.path.join(testfolder, name))
        if not os.path.exists(testsubf):
            os.makedirs(testsubf)

        # Read notebook in
        with open(tut, 'r') as fi:
            notebook = nbformat.read(fi, nbformat.NO_CONVERT)
            output, resources = pyexp.from_notebook_node(notebook, )

            # Write it out in .py format
            with open(os.path.join(testsubf, 'test_{}.py'.format(name)), 'w') as fo:
                fo.writelines(format_script(output, testsubf, name))

    return pytest.main(['-k', 'test_tutorial', testfolder])


def format_script(script, outdir, name):
    # 1. Wrap it all in a test function for pytest
    # 2. Indent everything one tab to be inside the function

    splitt = script.split('\n')
    lines = ['from htmd.ui import *\n',
             'def test_tutorial_{}():\n'.format(name.replace('-', '_')),
             '\tgetCurrentViewer(dispdev="text")\n',
             '\tnp.random.seed(0)\n',
             '\timport os\n',
             '\tos.chdir("{}")\n'.format(outdir)]
    for l in splitt:
        if not l.startswith('get_ipython(') and not l.startswith('from htmd.ui import *'):
            l = l.replace('webgl', 'vmd')
            l = l.replace('ngl', 'vmd')
            lines.append('\t' + l + '\n')
    return lines

