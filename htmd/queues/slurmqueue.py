from jobqueues.slurmqueue import SlurmQueue
import logging


logger = logging.getLogger(__name__)

logger.warning('Please do not import from htmd.queues.slurmqueue. This will be deprecated. '
                'Use jobqueues.slurmqueue instead.')