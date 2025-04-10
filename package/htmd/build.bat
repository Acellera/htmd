%PYTHON% -m build -v -w --no-isolation
for %%F in (dist\*.whl) do %PYTHON% -m pip install "%%F"
