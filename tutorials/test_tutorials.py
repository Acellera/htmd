import nbformat
from nbconvert import PythonExporter
from glob import glob
import os
import pytest
import htmd


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

    return pytest.main([testfolder])


def format_script(script, outdir, name):
    # 1. Need to remove the * imports. Python doesn't allow * imports in functions
    # 2. Wrap it all in a test function for pytest
    # 3. Indent everything one tab to be inside the function
    import htmd.ui
    htmdimports = [x for x in htmd.ui.__dict__.keys() if not x.startswith('_')]

    splitt = script.split('\n')
    lines = ['def test_{}():\n'.format(name.replace('-', '_')),
             '\tfrom htmd import {}\n'.format(', '.join(htmdimports)),
             '\timport matplotlib\n',
             '\tmatplotlib.use("Agg")\n',
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

