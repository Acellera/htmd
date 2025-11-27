# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
for d in ./*/; do
    cd $d
    name=`ls *.xtc`
    echo $name
    rm $name
    touch -d "1 year ago" $name
    cd ..
done
