This repository contains DSOs for Linux in htmd/lib/Linux
Source to these is in the restricted access repo Acellera/htmdlib, only available to Acellera-approved developers.
If you have access to that repo and want to build the DSOs, do the following:

cd [root of htmd checkout]
git clone ssh://git@github.com/acellera/htmdlib --depth 1
htmdlib/C/build.sh $PWD/htmd/lib/Linux/

