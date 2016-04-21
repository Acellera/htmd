# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import tempfile
import logging

logger = logging.getLogger(__name__)


def tempname(suffix='', create=False):
    if create:
        file = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    else:
        file = tempfile.NamedTemporaryFile(delete=True, suffix=suffix)
    file.close()
    return file.name





def _testDeprecation():
    issueDeprecationWarning("Please ignore this first warning")
    issueDeprecationWarning("This second warning should not appear")



_issuedDeprecationWarnings = {}


def issueDeprecationWarning(msg):
    """ Issue a deprecation warning, only once.

    Parameters
    ----------
    msg : str
        The message.

    """

    import inspect
    caller=inspect.stack()[1][3]
    if not caller in _issuedDeprecationWarnings:
        _issuedDeprecationWarnings[caller]=1
        logger.warning("Deprecation warning (%s). The function is obsolete and will be removed in a future version! Additional info: %s" % (caller,msg) )


if __name__ == "__main__":
    import doctest
    doctest.testmod()

    _testDeprecation()

