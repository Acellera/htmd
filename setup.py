# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import setuptools
import subprocess

try:
    version = (
        subprocess.check_output(["git", "describe", "--abbrev=0", "--tags"])
        .strip()
        .decode("utf-8")
    )
except Exception as e:
    print("Could not get version tag. Defaulting to version 0")
    version = "0"

with open("requirements.txt", "r") as f:
    requirements = f.read().splitlines()


if __name__ == "__main__":
    with open("README.md", "r") as fh:
        long_description = fh.read()

    setuptools.setup(
        name="acellera-htmd",
        version=version,
        description="High throughput molecular dynamics (HTMD)",
        long_description=long_description,
        zip_safe=False,
        url="https://www.htmd.org",
        maintainer="Acellera Ltd",
        maintainer_email="info@acellera.com",
        entry_points={
            "console_scripts": [
                "htmd_register   = htmdx.cli:htmd_do_register",
                "activate_license  = htmdx.cli:main_activate",
            ]
        },
        packages=setuptools.find_packages(include=["htmd*"], exclude=["data"]),
        package_data={"htmd": ["share/*/*/*/*", "logging.ini"]},
        install_requires=requirements,
    )
