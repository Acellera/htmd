from moleculekit.projections.metricdistance import MetricDistance, MetricSelfDistance
import logging


logger = logging.getLogger(__name__)

logger.warning('Please do not import from htmd.projections.metricdistance. This will be deprecated. '
                'Use moleculekit.projections.metricdistance instead.')