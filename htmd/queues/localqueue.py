from jobqueues.localqueue import LocalCPUQueue, LocalGPUQueue
import logging


logger = logging.getLogger(__name__)

logger.warning('Please do not import from htmd.queues.localqueue. This will be deprecated. '
                'Use jobqueues.localqueue instead.')