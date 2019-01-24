import os

from shutil import rmtree
from datetime import datetime

from dolfin import *
from .tissuebase import TissueBase
from .util import SmartFunction

from mpi4py import MPI

from json import JSONEncoder

class DolfinEncoder(JSONEncoder):
	def default(self, o):
		return str(o)

class Simulator(object):
	""" Handles running the simulation of the problem defined by a :class:`morphosolver.core.tissuebase.TissueBase` object
	
	Attributes:
		tissue (:class:`morphosolver.core.tissuebase.TissueBase`): tissue to simulate
		outputDir (:obj:`str`): path to an output directotry, defaults to the current workind directory 
	"""
	def __init__(self, tissue, outputDir = os.getcwd()):

		if MPI.COMM_WORLD.Get_rank() == 0:
			self.outputPath = os.path.join( outputDir, tissue.__class__.__name__, datetime.now().isoformat() )
			os.makedirs(self.outputPath)
		else:
			self.outputPath = None

		self.outputPath = MPI.COMM_WORLD.bcast(self.outputPath, root = 0)

		print("outputPath:", self.outputPath)

		try:
			self.tissue = tissue

			self.tissue.init()
			self.registerParams()

			paramFile = open( os.path.join(self.outputPath, "params.json"), "w" )
			paramFile.write( DolfinEncoder().encode(self.tissue.params) )
			paramFile.close()

		except Exception as e:
			rmtree(self.outputPath)
			raise e

	def run(self, dt, n, growthOn = True):
		""" Runs the simulation
		
		Args:
			dt (:obj:`float`): A time-step
			n (:obj:`int`): Number of time-steps
			growthOn (:obj:`bool`): Enable growth. Defaults to True
		"""
		print("Time {0:.2f}\n".format(self.tissue.t))
		print("\tSaving initial state...")
		self.saveParams()
		for i in range(n):
			
			print("Time {0:.2f}\n".format(self.tissue.t + dt))
			
			self.tissue.t += dt

			if growthOn:
				print("\tSolving growth...\n")
				self.tissue.solve(dt)

			print("\tUpdating state...\n")
			self.tissue.update(dt)

			print("\tSaving parameters...")
			self.saveParams()
			
	def saveParams(self):
		self.meshFile << (self.tissue.mesh, self.tissue.t)

		for paramName in self.tissue.params:
			if isinstance(self.tissue.params[paramName], SmartFunction):
				self.tissue.params[paramName].save(self.tissue.t)

		for attr in dir(self.tissue):
			if isinstance(getattr(self.tissue, attr), SmartFunction):
				getattr(self.tissue, attr).save(self.tissue.t)

	def registerParams(self):
		self.tissue.G = SmartFunction("G", self.tissue.W0, function = self.tissue.G)
		self.meshFile = File( os.path.join(self.outputPath, "mesh.pvd") )

		for paramName in self.tissue.params:		
			if isinstance(self.tissue.params[paramName], Function):
				if not isinstance(self.tissue.params[paramName], SmartFunction):
					function = self.tissue.params[paramName]
					if (function.value_rank() == 0):
						self.tissue.params[paramName] = SmartFunction(paramName, self.tissue.C0, function = function)
					elif (function.value_rank() == 1):
						self.tissue.params[paramName] = SmartFunction(paramName, self.tissue.V0, function = function)
					elif (function.value_rank() == 2):
						self.tissue.params[paramName] = SmartFunction(paramName, self.tissue.W0, function = function)
					else:
						raise(ValueError, "Invalid function.value_rank() {0}".format(function.value_rank()))
				
				self.tissue.params[paramName].setOutputPath(self.outputPath)

		for attr in dir(self.tissue):
			if isinstance( getattr(self.tissue, attr), SmartFunction ):
				getattr(self.tissue, attr).setOutputPath(self.outputPath)