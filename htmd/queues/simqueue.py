from jobqueues.simqueue import SimQueue, RetrieveError, SubmitError, InProgressError, ProjectNotExistError
import logging


logger = logging.getLogger(__name__)

logger.warning('Please do not import from htmd.queues.simqueue. This will be deprecated. '
                'Use jobqueues.simqueue instead.')
