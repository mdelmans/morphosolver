import os

from dolfin import *

from ...core.util import SmartFunction
from .gf import Gf

class AnisoGf(Gf):
	def __init__(self, tissue):
		super().__init__(tissue)

		self.a		= SmartFunction( "anisoGf_a",		functionSpace = self.tissue.C )
		self.g		= SmartFunction( "anisoGf_g",		functionSpace = self.tissue.C )
		self.gamma	= SmartFunction( "anisoGf_gamma",	functionSpace = self.tissue.V )

	def setOutputPath(self, outputPath):
		super().setOutputPath(outputPath)
		
		self.a.setOutputPath(outputPath)
		self.g.setOutputPath(outputPath)
		self.gamma.setOutputPath(outputPath)

	def projectAGG(self, a = Constant(0.0), g = Constant(0.0), gamma = Constant([0,0,1])):
		self.a.project(a)
		self.g.project(g)
		self.gamma.project( gamma / sqrt( inner(gamma, gamma) ) )

		g0 = g * ( 1 - a	)
		g1 = g * ( 1 + 2*a	)

		self.project( Identity(3) * g0 + (g1-g0) * outer( gamma, gamma ) )

	def save(self, t):
		self.file << (self, t)
		self.a.save(t)
		self.g.save(t)
		self.gamma.save(t)

	# def setGf(self, a = Constant(0.0), g = Constant(0.0), gamma = Constant([0,0,1]) ):
	# 	self.a.project(a)
	# 	self.g.project(g)
	# 	self.gamma.project( gamma / sqrt( inner(gamma, gamma) ) )

	# 	g0 = g * ( 1 - a	)
	# 	g1 = g * ( 1 + 2*a	)

	# 	self.gf.project( Identity(3) * g0 + (g1-g0) * outer( gamma, gamma ) )

# class Anisotropy:
# 	def setGrowth(self, a = Constant(0.0), g = Constant(0.0), gamma = Constant([0,0,1])):
# 		self.a = a
# 		self.g = g
# 		self.gamma = gamma

# 	def updateG(self, dt):
# 		g = project(self.g, self.C)
# 		a = project(self.a, self.C)
		
# 		gamma = project(self.gamma / sqrt( inner(self.gamma, self.gamma) ), self.V)

# 		g0 = g * ( 1 -		a )
# 		g1 = g * ( 1 + 2 *	a )
		
# 		gf = project( Identity(3) * g0 + (g1-g0) * outer( gamma, gamma ), self.W )

# 		self.G.project( (Function(self.W0, gf.vector()) * dt + Identity(3)) * self.G )