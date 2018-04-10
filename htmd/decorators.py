# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import functools
import inspect
import warnings


class _Deprecated(object):
    """
    Flags a class or function as deprecated

    Parameters
    ----------
    version : str
        The HTMD version at which the class or function was deprecated.
    newname : str
        The replacement function, with complete Namespace.

    Examples
    --------
    >>> @_Deprecated('0.0.0')
    ... class ExampleClass:
    ...     @_Deprecated('0.0.1')
    ...     def example_method(self):
    ...         pass
    ...     pass
    >>> ExampleClass.__doc__
    '\\n    .. warning:: Deprecated since version 0.0.0. \\n                '
    >>> ExampleClass().example_method.__doc__
    '\\n        .. warning:: Deprecated since version 0.0.1. \\n                '
    >>> ExampleClass() # doctest:+ELLIPSIS
    <__main__.ExampleClass object ...>
    >>> ExampleClass().example_method # doctest:+ELLIPSIS
    <bound method ExampleClass.example_method of <__main__.ExampleClass object ...>>
    """
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

if __name__ == "__main__":
    import doctest

    doctest.testmod()
