#!/bin/sh
# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
acemd --device 0 input > log.0 &
acemd --device 1 input > log.1 &
acemd --device 2 input > log.2 &
acemd --device 3 input > log.3 &
wait
