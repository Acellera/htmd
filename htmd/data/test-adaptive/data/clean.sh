for d in ./*/; do
    cd $d
    name=`ls *.xtc`
    echo $name
    rm $name
    touch -d "1 year ago" $name
    cd ..
done
