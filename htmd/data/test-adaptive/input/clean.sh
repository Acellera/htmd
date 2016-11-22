for d in ./*/; do
    cd $d
    rm ./*
    touch -d "1 year ago" structure.pdb
    cd ..
done
