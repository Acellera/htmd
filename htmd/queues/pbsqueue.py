from jobqueues.pbsqueue import PBSQueue
import logging


logger = logging.getLogger(__name__)

logger.warning('Please do not import from htmd.queues.pbsqueue. This will be deprecated. '
                'Use jobqueues.pbsqueue instead.')