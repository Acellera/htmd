# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

# This file is part of PyEMMA.
#
# Copyright (c) 2015, 2014 Computational Molecular Biology Group, Freie Universitaet Berlin (GER)
#
# PyEMMA is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

'''
Created on 24.04.2015

@author: marscher
'''

from __future__ import absolute_import
import sys
import time

__all__ = ('is_interactive_session', 'show_progressbar')


def __attached_to_ipy_notebook():
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            from IPython.html.widgets import IntProgress
            IntProgress(100)
        except:
            return False
        else:
            return True


def __is_interactive():
    # started by main function or interactive from python shell?
    import __main__ as main
    return not hasattr(main, '__file__')


def __is_tty_or_interactive_session():
    is_tty = sys.stdout.isatty()
    is_interactive = __is_interactive()
    result = is_tty or is_interactive
    return result

#ipython_notebook_session = __attached_to_ipy_notebook()
ipython_notebook_session = False
""" are we running an interactive IPython notebook session """

is_interactive_session = __is_tty_or_interactive_session()
""" do we have a tty or an interactive session? """

if ipython_notebook_session:
    from IPython.display import display
    from IPython.html.widgets import IntProgress, Box, Text


def hide_widget(widget):
    widget.close()


def hide_progressbar(bar):
    if ipython_notebook_session and hasattr(bar, 'widget'):
        from threading import Timer
        timeout = 2
        Timer(timeout, hide_widget, args=(bar.widget, )).start()
    else:
        sys.stdout.write("\n")


def show_progressbar(bar, show_eta=True, force=False):
    """ shows given bar either using an ipython widget, if in
    interactive session or simply use the string format of it and print it
    to stdout.

    Parameters
    ----------
    bar : instance of pyemma.util.progressbar.ProgressBar
    show_eta : bool (optional)

    """
    currtime = time.time()
    if bar.lastupdate is not None and not force:
        dtlastupdate = currtime - bar.lastupdate
    else:
        dtlastupdate = 9999
    #if not (str(config['show_progress_bars']) == 'True' and __is_tty_or_interactive_session() and dtlastupdate > 0.5):
    if not (__is_tty_or_interactive_session() and dtlastupdate > 0.5):
        return
    bar.lastupdate = currtime

    # note: this check ensures we have IPython.display and so on.
    if ipython_notebook_session:
        # create IPython widgets on first call
        if not hasattr(bar, 'widget'):
            box = Box()
            text = Text()
            progress_widget = IntProgress()

            box.children = [text, progress_widget]
            bar.widget = box
            widget = box

            # make it visible once
            display(box)

            # update css for a more compact view
            progress_widget._css = [
                ("div", "margin-top", "0px")
            ]
            progress_widget.height = "8px"
        else:
            widget = bar.widget

        # update widgets slider value and description text
        desc = bar.description
        desc += ':\t({}/{})'.format(bar.numerator, bar.denominator)
        if show_eta:
            desc += ':\tETA:' + bar._generate_eta(bar._eta.eta_seconds)
        assert isinstance(widget.children[0], Text)
        assert isinstance(widget.children[1], IntProgress)
        widget.children[0].placeholder = desc
        widget.children[1].value = bar.percent
    else:
        sys.stdout.write("\r" + str(bar))
        sys.stdout.flush()
