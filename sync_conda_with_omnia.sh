#!/usr/bin/env python3

import requests
import os
from subprocess import call

from binstar_client.utils import get_server_api
api=get_server_api()

os.chdir("/tmp")

try:
	os.mkdir("linux-64")
	os.mkdir("linux-32")
	os.mkdir("osx-64")
	os.mkdir("osx-32")
	os.mkdir("win-64")
	os.mkdir("win-32")
	os.mkdir("noarch")
except:
	pass

# Add packages to sync in this list here

for p in [ 
	"fftw3f",
	"openmm",
	"ambermini",
	"bhmm",
	"funcsigs",
	"mdtraj",
	"msmtools",
	"openbabel",
	"pint",
	"progress_reporter",
	"pyemma",
	"thermotools"
	]:
	omnia    = api.package("omnia", p );
	try:
		acellera = api.package("acellera", p );
	except:

		acellera = {"latest_version": "0" }

	if omnia["latest_version"] != acellera["latest_version"]:
		print("Syncing %s version %s.." %(p, omnia["latest_version"]) )
		for f in omnia["files"]:
			if f["version"] == omnia["latest_version"]:
				url = "https:" + f["download_url"]
				print("Downloading %s" %(url) )
				try:
					os.unlink("/tmp/package.bz2")
				except:
					pass
				call([ "curl", "-L", "-s",  url, "-o", f["basename"] ])
				print("Uploading.."  )
				try:
					os.getenv("ANACONDA_TOKEN_BASIC") 
					call([ "anaconda", "upload", "-t", os.getenv("ANACONDA_TOKEN_BASIC"),"-u",  "acellera", f["basename"] ])
				except:
					call([ "anaconda", "upload", "-u",  "acellera", f["basename"] ])
	else:
		print("Package %s up to date at version %s" % ( p, omnia["latest_version"] ) )
