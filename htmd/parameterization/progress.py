# (c) 2015-2019 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

class Progress:

    def __init__(self, job, num_stages):

        self._job = job
        self._stage = 0
        self._num_stages = num_stages

    def __call__(self, message):

        if self._job is not None:
            self._stage += 1
            self._job.setProgress(message, f'{self._stage}/{self._num_stages}')