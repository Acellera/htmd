# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
for d in ./*/; do
    cd $d
    rm ./*
    touch -d "1 year ago" structure.pdb
    cd ..
done
