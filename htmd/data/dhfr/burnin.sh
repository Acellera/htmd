#!/bin/sh
acemd --device 0 input > log.0 &
acemd --device 1 input > log.1 &
acemd --device 2 input > log.2 &
acemd --device 3 input > log.3 &
wait
