from setuptools import setup
import subprocess

version=subprocess.Popen(["git", "describe", "--tags"], stdout=subprocess.PIPE ).stdout.read().decode("utf8")
version=version.split("-")
version=version[0]

f=open( "package/htmd-deps/DEPENDENCIES", "r" )
deps=[ 'pyEMMA<2.3' ]
#for a in f.readlines():
#	deps.append(a.strip().split("=")[0]) #.decode("utf8"))

print("Version [%s]" % (version) )
print("Dependencies:" )
print(deps)

setup(name='htmd',
	version='0.15.11',
	description='HTMD',
	packages=['htmd', 'htmdx', 'htmd.molecule', 'htmd.parameterization', 'htmd.adaptive', 'htmd.apps', 'htmd.builder', 'htmd.clustering', 'htmd.progress', 'htmd.projections', 'htmd.protocols', 'htmd.qm', 'htmd.queues', 'htmd.data' ],
	install_requires=deps,
	zip_safe=False,
	url="https://www.htmd.org",
	maintainer="Acellera Ltd",
	maintainer_email="noreply@acellera.com",
	entry_points = {
		"console_scripts": [
			'htmdnb    = htmdx.cli:main_htmd_notebook',
			'htmd      = htmdx.cli:main_htmd',
			'htmd_register   = htmdx.cli:main_do_nothing',
			'parameterize = htmd.parameterization.cli:main_parameterize',
			'activate_license  = htmdx.cli:main_activate'
		]
	},
	include_package_data=True
	)

