import traceback
import versioneer
import sys

dep_file = sys.argv[1]

try:
    __version__ = versioneer.get_version()
except Exception:
    print(traceback.format_exc())
    print("Could not get version. Defaulting to version 0")
    __version__ = "0"

with open(dep_file, "r") as f:
    deps = [line for line in f if not line.startswith("#") and len(line.strip())]

# Fix conda meta.yaml
with open("package/htmd/meta.yaml", "r") as f:
    text = f.read()

text = text.replace("BUILD_VERSION_PLACEHOLDER", __version__)

text = text.replace(
    "DEPENDENCY_PLACEHOLDER",
    "".join(["    - {}\n".format(dep.strip()) for dep in deps]),
)

with open("package/htmd/meta.yaml", "w") as f:
    f.write(text)
