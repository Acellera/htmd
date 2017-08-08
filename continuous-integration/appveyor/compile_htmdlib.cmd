cd ..
git clone https://%GITHUB_HTMDLIB_USERNAME%:%GITHUB_HTMDLIB_PASSWORD%@github.com/Acellera/htmdlib --depth 1
choco install mingw
refreshenv
SET PATH=%PYTHON%;%PYTHON%\\Scripts;%PATH%
conda create -q -n py27 python=2.7
activate py27
conda install -c anaconda scons
cd htmdlib/C
scons --prefix=../../htmd/htmd/lib/Windows/
cd ../../htmd/
deactivate py27