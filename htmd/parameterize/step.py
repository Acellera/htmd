import inspect
import os
import logging
import shutil
import sys
import subprocess

log = logging.getLogger(__name__)



class StepSet:

	def __init__( self ):
		self._steps=[]
		self._steps_by_index={}


	def run( self, rootdir, param=None ):

		for i in self._steps:
			log.debug( "Step : [" + i.name + "]")
			log.debug( "-- Input files--")
			for j in i.inputlist:
				log.debug( j )
			log.debug( "-- Output files--")
			for j in i.outputlist:
				log.debug( j )

		log.info( "run(): Starting run in [" + rootdir + "]")
		for i in self._steps:
			directory = os.path.join( rootdir , os.path.basename(i.directory_name()) ) #"%03d-%s" % ( i.index, i.name ) )
			log.info( "run(): Step [" + i.name +"] in dir [" + directory +"]")
			if not os.path.exists( directory ):
				os.mkdir( directory )
			if not i.completed( directory ):
				self.stagein( i.inputlist, directory )
				log.info( "run(): Running step [" + i.name +"]" )
				errfile = os.path.join( directory, "error.txt" )
				if os.path.exists(errfile):
					os.unlink(errfile)

				print( "Running step [%03d] [%s]" % (i.index, i.name) )
				log.info( " Running step [" + str(i.index) + "] \t[" + i.name +"] " )

				sys.stdout.flush()
				i.run( directory, param=param )
				if os.path.exists(errfile):
					f=open( errfile, "r" )
					err=f.readlines()
					f.close()
					ret = "run(): Step yielded error: \n"
					for i in err:
						ret = ret + "\t" + i
					raise ValueError( ret )

				if not i.completed( directory ):
					raise ValueError( "run(): Step [" +  i.name + "] did not produce the expected output: " + str( i.missing(directory)) + " in directory " + directory )


	def add( self, obj ):
		for i in range(len(self._steps) ):
			if self._steps[i].index > obj.index:
				self._steps.insert(i, obj )
				return
		self._steps.append( obj )
		self._steps_by_index[ obj.index ]=obj
		obj.set = self

		for i in self._steps:
			for f in i.inputlist:
				t=f.split( ":" )
				idx=int(t[0])
				try:
					s = self._steps_by_index[idx]
				except:
					raise ValueError( "Step [" + str(idx) + "] is not defined" )

				if not f in s.outputlist:
					s.outputlist.append( t[1] )
					#log.debug("add(): Adding expected output file [" + t[1] + "] to step [" + s.name + "]" )

	def stagein( self, inputlist, directory ):
		import fnmatch

		rootdir = os.path.dirname( directory )
		for f in inputlist:
			t=f.split( ":" )
			filename=t[0]
			sourcestep=self._steps_by_index[int(t[0])]
			sourcedir = os.path.join( rootdir, sourcestep.directory_name() )
			sourcefile= os.path.join( sourcedir, t[1] )
			found = False
			for ff in os.listdir( sourcedir ):
				destname=t[1]
				if len(t)>2:
					destname=t[2]
				if fnmatch.fnmatch( ff, t[1]):  # Deal with matching globs
					log.debug("stagein(): Copying [" + sourcefile + "] to [" + destname + "]" )
					sf = os.path.join( sourcedir, ff )
					if len(t)>2:
						df = os.path.join( directory, destname )
					else:
						df = os.path.join( directory, ff )
					shutil.copyfile( sf, df )
					found = True

			if not found:
				raise ValueError( "stagein(): Expected output file [" + sourcefile + "] not found" )


class Step:
	def __init__( self, index, name, inputlist, outputlist  ):
		self.index     = index
		self.name      = name
		self.inputlist = inputlist
		if outputlist:
			self.outputlist= outputlist
		else:
			self.outputlist=[]

		self.set = None

	def run( self, directory, param=None ):
		log.info("Step [%s] doing nothing" % ( self.name ) )

	def directory_name(self):
		return "%03d-%s" % ( self.index, self.name )



	def missing( self, directory):
		import fnmatch
		import os
		missing=[];
		for i in self.outputlist:
			found = False
			for ff in os.listdir( directory ):
				if fnmatch.fnmatch(  ff, i ):
					found=True

			if not found:
				if not i in missing:
					missing.append(i);

		return missing



	def completed( self, directory ):
		import fnmatch
		import os
		for i in self.outputlist:
			found = False
			for ff in os.listdir( directory ):
				if fnmatch.fnmatch(  ff, i ):
					found=True

			if not found:
				log.info( "completed(): Expected outputfile [" + i +"] not present")
				return False
		return True


class ScriptStep(Step):
	def __init__( self, index, prefix, name, inputfiles, outputfiles=None ):
		root = inspect.getfile( ScriptStep )
		path = os.path.dirname(root) 
		path = os.path.join( path, "share" )
		path = os.path.join( path, "scripts" )
		path = os.path.join( path, prefix )
		path = os.path.join( path, "%03d-%s" % ( index, name ) )

		if not os.path.exists( path ) and not os.access( path, os.R_OK ):
			raise ValueError( "Script [" + path + "] not found"  )

		self._script = path

		super(ScriptStep, self).__init__( index, name, inputfiles, outputfiles  )

#	def dostep( self, directory ):
#		print(" Step [%03d-%s] doing nothing in [%s]" % ( self.index, self.name, directory )

	def run( self, directory, param=None ):
		import shutil
		import subprocess
		import stat
		dest = os.path.join( directory, "run.sh" )
		shutil.copyfile( self._script, dest)
		os.chmod( dest, stat.S_IEXEC | stat.S_IREAD | stat.S_IWRITE )
		pwd = os.getcwd()
		os.chdir( directory )
		ret=""
		err=None
		try:
			ret =subprocess.check_output( [ dest ], stderr = subprocess.STDOUT, shell=False, stdin = None )

		except subprocess.CalledProcessError as e:

			err=""
			err = err + "\n Failure in step   : " + self.name + "\n"
			err = err + "   Directory       : " + directory + "\n"
			err = err + "   Return value    : " + str(e.returncode) + "\n"
			err = err + "   Output          : " + "\n"
			err = err +  " -- " + "\n"
			for i in  e.output.decode("latin1").split( "\n" ):
				err = err +  i + "\n"
			err = err + " -- " + "\n"
			# return value was non-zero
		os.chdir( pwd )
		if err is not None:
			raise ValueError( err )


class QMScriptStep(Step):
	def __init__( self, index, prefix, name, inputfiles, outputfiles=None ):
		root = inspect.getfile( ScriptStep )
		path = os.path.dirname(root) 
		path = os.path.join( path, "share" )
		path = os.path.join( path, "scripts" )
		path = os.path.join( path, prefix )
		path = os.path.join( path, "%03d-%s" % ( index, name ) )

		if not os.path.exists( path ) and not os.access( path, os.R_OK ):
			raise ValueError( "Script [" + path + "] not found"  )

		self._script = path

		super(QMScriptStep, self).__init__( index, name, inputfiles, outputfiles  )

#	def dostep( self, directory ):
#		print(" Step [%03d-%s] doing nothing in [%s]" % ( self.index, self.name, directory )

	def _exec( self, cmd ):
		try:
			ret =subprocess.check_output( cmd , stderr = subprocess.STDOUT, shell=False, stdin = None )
		except subprocess.CalledProcessError as e:

			err=""
			err = err + "\n Failure in step   : " + self.name + "\n"
#			err = err + "   Directory       : " + directory + "\n"
			err = err + "   Return value    : " + str(e.returncode) + "\n"
			err = err + "   Output          : " + "\n"
			err = err +  " -- " + "\n"
			for i in  e.output.decode("latin1").split( "\n" ):
				err = err +  i + "\n"
			err = err + " -- " + "\n"
			# return value was non-zero
			if err is not None:
				raise ValueError( err )


	def run( self, directory, param=None ):
		import shutil
		import subprocess
		import stat
		import htmd.parameterize
		dest = os.path.join( directory, "run.sh" )
		shutil.copyfile( self._script, dest)
		os.chmod( dest, stat.S_IEXEC | stat.S_IREAD | stat.S_IWRITE )
		pwd = os.getcwd()
		os.chdir( directory )
		ret=""
		err=None

		self._exec( [ dest, "--prepare" ] )
		param.run_qm_jobs( directory )			
		self._exec( [ dest, "--complete"] )

		os.chdir( pwd )
		

