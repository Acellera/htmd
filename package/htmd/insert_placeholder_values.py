import yaml
import sys
import setuptools_scm

try:
    __version__ = setuptools_scm.get_version()
except Exception:
    print("Could not get version. Defaulting to version 0")
    __version__ = "0"

dep_file = sys.argv[1]

with open(dep_file, "r") as f:
    deps = [
        line.strip() for line in f if not line.startswith("#") and len(line.strip()) > 0
    ]

for i, dep in enumerate(deps):
    parts = [part.strip() for part in dep.split("#")]
    if len(parts) > 1 and parts[1].startswith("["):
        deps[i] = {"if": parts[1][1:-1], "then": parts[0]}

# Fix conda meta.yaml
with open("package/htmd/recipe_template.yaml", "r") as f:
    recipe = yaml.load(f, Loader=yaml.FullLoader)

recipe["package"]["version"] = __version__
recipe["requirements"]["run"] = deps

with open("package/htmd/recipe.yaml", "w") as f:
    yaml.dump(recipe, f)
