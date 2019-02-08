import sys
import os

from math import exp

from dolfin import *
from mshr import *

import numpy as np

from random import uniform

import matplotlib.pyplot as plt
import matplotlib.tri as tri

from morphosolver import TissueBase
from morphosolver import Simulator
from morphosolver import SmartFunction, Chemical, GenericDomain

from fenicstools.Interpolation import interpolate_nonmatching_mesh

# The function defying the state of the tissue, and the extent to which it is a source of inhibitor.
class StateFunction(UserExpression):
	def __init__(self, notches, k=11, **kwargs):
		self._notches = notches
		self._k = k
		
		super().__init__(**kwargs)
	def eval(self, value, x):
		closest = None
		minR = 1e30

		value[0] = 0.0

		if len(self._notches) > 0:

			for notch in self._notches:
				r = Point(x[0], x[1], 0).distance(notch)
				
				if r < minR:
					minR = r
					closest = notch

			centroid = sum(self._notches, Point(0,0)) / len(self._notches)
			dist = 2*centroid.distance(closest)
			
			value[0] = exp(-self._k * (minR/dist)**2 )

	def value_shape(self):
		return (1,)

	def update(self, tissue, **params):
		tissue.ut.set_allow_extrapolation(True)

		self._notches = []

		for notch in tissue.notches0:
			self._notches.append( notch + Point( tissue.ut( notch ) ) )

class InhibitorChemistry(object):
	def init(self, tissue, **params):
		self.params = params

		ci0 = self.params['ci0']
		co0 = self.params['co0']
		D   = self.params['D']

		tissue.ci = Chemical(tissue, "ci")
		tissue.co = Chemical(tissue, "co")
		tissue.co.setDiffusion( D )

		tissue.ci.project( ci0 )
		tissue.co.project( co0 )

	def update(self, tissue, dt):
		sf = tissue.stateField
		
		ci = tissue.ci.c0
		co = tissue.co.c0

		p	= self.params['p']
		ke	= self.params['ke']
		ki	= self.params['ki']
		d	= self.params['d']

		g  = tissue.gf.g

		tissue.ci.setReaction( sf * p  -  sf * ke * ci  +  (1-sf) * ki * co  -  (1-sf) * d * ci  - 3*g*ci	)
		tissue.co.setReaction(            sf * ke * ci  -  (1-sf) * ki * co        				 - 3*g*co   )

		tissue.ci.update(dt)
		tissue.co.update(dt)

# A Mix-in class for generating gemma mesh
class GemmaMeshGenerator(object):
	def init(self, tissue, **params):
		
		a = params['a']
		b = params['b']
		sharpness = params['sharpness']
		num = params['num']
		refinement = params['refinement']
		nd = params['nd']
		fileName = params['fileName']

		theta1 = np.linspace(2*np.pi/num, np.pi*(1 - 2/num), num = num)
		theta2 = np.linspace(np.pi*(1 + 2 / num), np.pi*(2 - 2/num), num = num)
		
		r1 = a*np.power( np.cos( -theta1 + 0.5*np.pi  ), 0.25 )  + b
		r2 = a*np.power( np.cos(  theta2 + 0.5*np.pi  ), 0.25 )  + b
		
		x1 = r1*np.cos(theta1)
		x2 = r2*np.cos(theta2)

		y1 = r1*np.sin(theta1)
		y2 = r2*np.sin(theta2)

		points = [ Point(x,y) for x,y in zip( x1, y1 ) ] + [ Point(x,y) for x,y in zip(x2, y2) ]

		polygon = Polygon(points)
		extruded = Extrude2D(polygon, 0.1)

		tissue.notches0 = [Point( x1[0] - nd,  0), Point(-x1[0] + nd,  0)]

		if fileName and os.path.isfile(fileName):
			return Mesh(fileName)
		
		if fileName:
			mesh = generate_mesh(extruded, refinement)
			meshFile = File(fileName)
			meshFile << mesh
		
		return mesh

class SliceMeshGenerator(object):
	def init(self, tissue, **params):
		tissue.notches0 = params['notches']

		fileName = params['fileName']

		if fileName and os.path.isfile(fileName):
			return Mesh(fileName)

		else:
			originalMesh = params['mesh']
			width = params['width']

			domains = MeshFunction('size_t', originalMesh, originalMesh.topology().dim())
			domains.set_all(0)

			sliceDomain = GenericDomain(lambda x, onBoundary: abs(x[0]) <= 0.5 * width )
			sliceDomain.mark(domains, 1)

			mesh = SubMesh(originalMesh, domains, 1)

			if fileName:
				meshFile = File(fileName)
				meshFile << mesh

			return mesh 

class NotchGrowth(TissueBase):
	def __init__(self, meshGenerator, chemistry, params):
		params['C1'] = 100

		
		mesh =  meshGenerator.init(self, **params['meshParams'])

		self.chemistry = chemistry

		super().__init__(mesh, params)
		self.stateField = StateFunction(self.notches0, element=self.C.ufl_element())

	def init(self):
		self.gf = AnisoGf(self)
		self.sf = SmartFunction("sf", function = project(self.stateField, self.C) )

		self.chemistry.init(self, **self.params['chemistryParams'])

	def update(self, dt):
		self.stateField.update(self)
		self.chemistry.update(self, dt)

if __name__ == "__main__":
	
	# Run diffusion
	meshParams = dict(a = 0.5, b = 0.7,  sharpness = 0.5, num = 200, refinement = 50, nd = 0.1, fileName = 'mesh.xml')
	chemistryParams = dict(ci0 = Constant(0.0), co0 = Constant(0.0), D = 1.0, p = 0.05, ke = 0.5, ki = 0.5, d = 0.05, sk = 0.05)
	
	params = dict(meshParams = meshParams, chemistryParams = chemistryParams)

	chemistry = InhibitorChemistry()
	meshGenerator = GemmaMeshGenerator()

	model0 = NotchGrowth(meshGenerator, chemistry, params = params)
	simulator0 = Simulator(model0, outputDir = sys.path[0])
	simulator0.run( float(sys.argv[1]), int(sys.argv[2]), growthOn = False )

	# Cut and wait
	cutMeshParams = dict( mesh = model0.mesh, width = 0.8, notches = [], fileName = 'cutMesh.xml')
	params = dict(meshParams = cutMeshParams, chemistryParams = chemistryParams)
	meshGenerator = SliceMeshGenerator()

	model1 = NotchGrowth(meshGenerator, chemistry, params = params)
	simulator1 = Simulator(model1, outputDir = sys.path[0])

	model1.ci.project( interpolate_nonmatching_mesh( model0.ci, model1.C) )
	model1.co.project( interpolate_nonmatching_mesh( model0.co, model1.C) )

	simulator1.run( float(sys.argv[1]), int(sys.argv[3]), growthOn = False )

	# Establish new notches
	model1.notches0 = [ Point(0,0.4,0), Point(0,-0.4,0) ]
	model1.stateField.notches = [ Point(0,0.4,0), Point(0,-0.4,0) ]
	simulator1.run( float(sys.argv[1]), int(sys.argv[4]), growthOn = False )