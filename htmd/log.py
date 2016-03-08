# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import logging


class HTMDlogger:
    logger = None

    def __init__(self, name):
        if HTMDlogger.logger is None:
            formatter = logging.Formatter(fmt='%(asctime)s - %(levelname)s - %(module)s - %(message)s')
            handler = logging.StreamHandler()
            handler.setFormatter(formatter)

            HTMDlogger.logger = logging.getLogger(name)
            HTMDlogger.logger.setLevel(logging.DEBUG)
            HTMDlogger.logger.addHandler(handler)

    def __getattr__(self, name):
        from IPython.core.debugger import Tracer
        Tracer()()
        return getattr(self.logger, name)

