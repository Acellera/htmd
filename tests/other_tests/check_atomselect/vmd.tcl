mol new minimal.pdb  filebonds off autobonds off
set a [ atomselect "top" "segname P1 and not protein" ]
$a num

mol new minimal.pdb  filebonds on autobonds off
set a [ atomselect "top" "segname P1 and not protein" ]
$a num


mol new minimal.pdb  filebonds off autobonds on
set a [ atomselect "top" "segname P1 and not protein" ]
$a num

mol new minimal.pdb  filebonds on autobonds on
set a [ atomselect "top" "segname P1 and not protein" ]
$a num
quit
