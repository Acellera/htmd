# (c) 2015-2026 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import subprocess
import sys


def _safe_run(cmd):
    try:
        return subprocess.check_output(cmd, shell=True).decode("utf8").strip()
    except subprocess.CalledProcessError:
        return ""


if sys.argv[1] == "tag":
    output = _safe_run("git describe --tags")
    tag = output.split("-")[0] if output else ""
    print(tag)

if sys.argv[1] == "branch":
    output = _safe_run("git rev-parse --abbrev-ref HEAD")
    if output.startswith("master"):
        print("latest")
    elif output.startswith("rel-"):
        print("stable")
    else:
        print(output or "latest")
