from dolfin import *
from ...core.util import SmartFunction


class Gf(SmartFunction):
	def __init__(self, tissue):
		self.tissue = tissue
		super().__init__("Gf", functionSpace = self.tissue.W, functionSpace0 = self.tissue.W0)

	def updateG(self, dt):
		self.tissue.G.project( (self.initial * dt + Identity(3)) * self.tissue.G )

# class Gf(object):
# 	def __init__(self, tissue):
# 		self.tissue = tissue
# 		self.gf = SmartFunction("Gf", functionSpace = self.tissue.W)

# 	def setGf(self, gf):
# 		self.gf.project(gf)

# 	def updateG(self, dt):
# 		self.tissue.G.project( (Function(self.tissue.W0, self.gf.vector()) * dt + Identity(3)) * self.tissue.G )