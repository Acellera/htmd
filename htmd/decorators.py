# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
import functools
import inspect
import warnings


class _Deprecated(object):
    def __init__(self, version, newname=None):
        self.version = version
        self.newname = newname

    def __call__(self, wrapped):

        if inspect.isfunction(wrapped):
            objtype = 'function/method'
            docprefix = \
                """
        .. warning:: Deprecated since version {}. {}
                """.format(self.version, 'Use :func:`{}` instead.'.format(self.newname) if self.newname else '')
        elif inspect.isclass(wrapped):
            objtype = 'class'
            docprefix = \
                """
    .. warning:: Deprecated since version {}. {}
                """.format(self.version, 'Use :func:`{}` instead.'.format(self.newname) if self.newname else '')
        else:
            raise TypeError(type(wrapped))

        message = "The `{}` {} has been deprecated " \
                  "since version {}. {}".format(wrapped.__name__, objtype,
                                                self.version,
                                                'Use `{}` instead.'.format(self.newname) if self.newname else '')

        wrappingclass = isinstance(wrapped, type)
        if wrappingclass:
            wrappedfunc = wrapped.__init__
        else:
            wrappedfunc = wrapped

        wrapped.__doc__ = docprefix + (wrapped.__doc__ or '')

        @functools.wraps(wrappedfunc)
        def new_func(*args, **kwargs):
            warnings.warn(message, category=DeprecationWarning, stacklevel=2)
            return wrappedfunc(*args, **kwargs)

        if wrappingclass:
            wrapped.__init__ = new_func
            return wrapped
        else:
            return new_func
