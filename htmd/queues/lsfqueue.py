from jobqueues.lsfqueue import LsfQueue
import logging


logger = logging.getLogger(__name__)

logger.warning('Please do not import from htmd.queues.lsfqueue. This will be deprecated. '
                'Use jobqueues.lsfqueue instead.')