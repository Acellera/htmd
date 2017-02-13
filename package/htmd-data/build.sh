# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
printenv

if [ "$?" != "0" ]; then
	echo "Error: Build failed"
	exit 1
fi


rm -rf $(find htmd -type d -name __pycache__)

	DIR="$SP_DIR"
	echo "Installing into $DIR"


	if [ "$DIR" != "" ]; then
		mkdir -p "$DIR"
	fi
	if [ -e "$DIR" ]; then
		pwd
		ls htmd/data
    mkdir -p $DIR/htmd
		cp -r htmd/data  $DIR/htmd/
    rm -rf $(find "$DIR/htmd" -name .git -type d) 

 	else
		echo "Error: SP_DIR not defined"
		exit 1
	fi

chmod -R a+rX .
chmod -R a+rX .



exit 0
