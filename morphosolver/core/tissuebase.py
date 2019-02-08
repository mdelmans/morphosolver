import os

# from dolfin import *
from dolfin import UnitCubeMesh, Mesh
from dolfin import Constant, FunctionSpace, VectorFunctionSpace, TensorFunctionSpace
from dolfin import project, inner, derivative, ALE, Identity
from dolfin import TestFunction, TrialFunction, Function
from dolfin import FacetNormal
from dolfin import SpatialCoordinate
from dolfin import parameters
from dolfin import dx, ds
from dolfin import NonlinearVariationalProblem, NonlinearVariationalSolver

import mshr

from ufl import inv, nabla_grad, min_value, grad, sym
import matplotlib.pyplot as plt

from .util import GenericDomain, SmartFunction

class TissueBase(object):
	"""
	TissueBase defines the problem geometry and mechanics. It is a base class, where :func:`~init` and :func:`~update`\
	functions should be overridden by the derived class.

	Args:
		mesh0 (:class:`dolfin.Mesh`): A mesh of the initial condition
		params (:obj:`dict`): A dictionary of parameters
	"""

	def __init__(self, mesh0 = UnitCubeMesh(8,8,8), params={}):
		parameters['form_compiler']['representation'] = 'uflacs'
		parameters['form_compiler']['optimize'] = True
		parameters['form_compiler']['quadrature_degree'] = 4

		self.mesh0 = Mesh(mesh0)
		self.mesh = Mesh(mesh0)

		if not 'C1' in params:
			params['C1'] = 100

		self.params = params

		self.b = Constant((0.0, 0.0, 0.0))
		self.h = Constant((0.0, 0.0, 0.0))

		self.C0 = FunctionSpace			(self.mesh0, "Lagrange", 2)
		self.V0 = VectorFunctionSpace	(self.mesh0, "Lagrange", 1)
		self.W0 = TensorFunctionSpace	(self.mesh0, "Lagrange", 1)

		self.C = FunctionSpace			(self.mesh, "Lagrange", 2)
		self.V = VectorFunctionSpace	(self.mesh, "Lagrange", 1)
		self.W = TensorFunctionSpace	(self.mesh, "Lagrange", 1)
		
		self.G = project( Identity(3), self.W0 )

		self.ut = Function(self.V0)

		self.du = TestFunction(self.V0)
		self.w 	= TrialFunction(self.V0)

		self.n0 = FacetNormal(self.mesh0)

		self.v = Function(self.V)

		self.t = 0.0


	def init(self):
		"""
		Virtual initialisation function. It ss only run once before the start of the simulation.
		"""
		pass

	def update(self,dt):
		"""
		Virtual iupdate function. Is run at every time-step of the simulation before solving the growth mechanics\
		by the call to the :func:`~solve` function.

		Args:
			dt (:obj:`float`): A time-step, passed by the :class:`~morphosolver.core.simulator.Simulator` object.
		"""
		pass

	def eps(self, u):
		return sym( nabla_grad(u) )

	def F(self):
		return Identity(3) + grad(self.u + self.ut)

	def A(self):
		return self.F() * inv(self.G)

	def S(self):
		C1 = self.params['C1']
		return C1*( self.A() - Identity(3) )

	def solve(self, dt):	
		self.u	= Function		(self.V0)
		self.w	= TestFunction	(self.V0)
		self.du	= TrialFunction	(self.V0)

		x = SpatialCoordinate(self.mesh0)

		L = inner( self.S(), self.eps(self.w) )*dx(degree=4)\
		- inner( self.b, self.w )*dx(degree=4)\
		- inner( self.h, self.w )*ds(degree=4)\
		+ inner( 1e-6*self.u, self.w )*ds(degree=4)\
		- inner( min_value(x[2]+self.ut[2]+self.u[2], 0) * Constant((0,0,-1.0)), self.w )*ds(degree=4)

		a = derivative(L, self.u, self.du)

		problem = NonlinearVariationalProblem(L, self.u, bcs=[], J=a)
		solver = NonlinearVariationalSolver(problem)

		solver.solve()

		self.ut.vector()[:] = self.ut.vector()[:] + self.u.vector()[:]

		ALE.move(self.mesh, Function(self.V, self.u.vector()) )

		self.v.vector()[:] = self.u.vector()[:] / dt
		self.n = FacetNormal(self.mesh)

	def plotMesh(self):
		"""
		Virtual iupdate function. Is run at every time-step of the simulation before solving the growth mechanics\
		by the call to the :func:`~solve` function.
		"""
		plt.figure()
		plot(self.mesh)