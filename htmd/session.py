# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

#import dill
import pickle
import logging
logger = logging.getLogger(__name__)


def htmdsave(fname, varnames=None):
    """ Saves the workspace to a file

    Parameters
    ----------
    fname : str
        The file in which to save the workspace
    varnames : list
        A list of variables to save. Default: saves all variables.
    """
    from types import ModuleType, FunctionType
    import __main__ as _main_module
    import inspect
    mydict = dict()
    maindict = _main_module.__dict__.copy()
    bannedvars = ['In', 'Out', 'exit', 'quit', 'get_ipython']

    if varnames is None:
        logger.warning('Save will ignore all variables starting with _ and following variable names {}'.format(bannedvars))
    logger.warning('Save will duplicate data which is stored in various pointers or variables. Take care!')

    if varnames is None:
        for k in maindict:
            if k[0] != '_' and k not in bannedvars and \
                    not isinstance(maindict[k], ModuleType) and \
                    not isinstance(maindict[k], FunctionType) and \
                    not inspect.isclass(maindict[k]):
                try:
                    pickle.dumps(maindict[k])  # This probably slows down everything
                    mydict[k] = maindict[k]
                except:
                    logger.debug('Skipping variable: {}'.format(k))
    else:
        if not isinstance(varnames, list):
            varnames = [varnames]
        for k in varnames:
            mydict[k] = maindict[k]

    logger.debug('Saving variables: {}'.format(list(mydict.keys())))
    f = open(fname, 'wb')
    pickle.dump(mydict, f)
    f.close()


def htmdload(fname):
    """ Loads a previously saved workspace from a file.

    Parameters
    ----------
    fname : str
        The location of the file containing the saved workspace.
    """
    import __main__ as _main_module

    f = open(fname, 'rb')
    mydict = pickle.load(f)
    f.close()

    # Inserting all variables into the main dictionary
    _main_module.__dict__.update(mydict)


'''def save(fname):
    from htmd.metricdata import MetricData
    from htmd.model import Model
    from htmd.kinetics import Kinetics

    import __main__ as _main_module
    myobj = dict()
    for k in _main_module.__dict__:
        if isinstance(_main_module.__dict__[k], MetricData) or isinstance(_main_module.__dict__[k], Model):
            myobj[k] = _main_module.__dict__[k]
        # TODO: Support lists and np.arrays of class objects

    f = open(fname, 'wb')
    pickle.dump(myobj, f)
    f.close()
    #dill.dump_session(filename=filename)
'''

'''
def save(dir)
    from htmd.metricdata import MetricData
    from htmd.model import Model
    from htmd.kinetics import Kinetics

    if not path.isdir(dir):
        os.makedirs(dir)

    import __main__ as _main_module
    myobj = dict()
    for k in _main_module.__dict__:
        if isinstance(_main_module.__dict__[k], MetricData):
            _main_module.__dict__[k].save(path.join(dir, 'MetricData_{}.dat'.format(k)))
        elif isinstance(_main_module.__dict__[k], Model):
            _main_module.__dict__[k].save(path.join(dir, 'Model_{}.dat'.format(k)))

def load(dir):
    from htmd.metricdata import MetricData
    from htmd.model import Model
    from htmd.kinetics import Kinetics
    from glob import glob

    if not path.isdir(dir):
        raise NameError('Directory {} does not exist'.format(dir))

    import __main__ as _main_module

    regex = re.compile('(\w+?)_(\w+).dat')

    datfiles = glob(path.join(dir, '*.dat'))
    for d in datfiles:
        res = regex.search(path.basename(d))
        if not res:  # If we are running on top of adaptive, use the first name part for the next sim name
            continue

        classname = res.group(1)
        varname = res.group(2)
        print(classname, varname)

        if classname == 'MetricData':
            _main_module.__dict__[varname] = MetricData()
        elif classname == 'Model':
            _main_module.__dict__[varname] = Model()
        _main_module.__dict__[varname].load(d)'''
