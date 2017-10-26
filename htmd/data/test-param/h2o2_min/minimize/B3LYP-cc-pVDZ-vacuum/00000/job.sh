#!/bin/bash


trap "touch /tmp/stef/h2o2_min/minimize/B3LYP-cc-pVDZ-vacuum/00000/htmd.queues.done" EXIT SIGTERM


trap "touch /tmp/stef/h2o2_min/minimize/B3LYP-cc-pVDZ-vacuum/00000/htmd.queues.done" EXIT SIGTERM

cd /tmp/stef/h2o2_min/minimize/B3LYP-cc-pVDZ-vacuum/00000
/tmp/stef/h2o2_min/minimize/B3LYP-cc-pVDZ-vacuum/00000/run.sh