import os

from dolfin import UnitCubeMesh, Mesh
from dolfin import Constant, FunctionSpace, VectorFunctionSpace, TensorFunctionSpace
from dolfin import project, inner, derivative, ALE
from dolfin import TestFunction, TrialFunction, Function
from dolfin import FacetNormal
from dolfin import SpatialCoordinate

import mshr

from ufl import inv, nabla_grad, min_value
import matplotlib.pyplot as plt

from .util import GenericDomain, SmartFunction

class TissueBase(object):
	"""
	TissueBase defines the problem geometry and mechanics. It is a base class, where :func:`~init` and :func:`~update`\
	functions should be overridden by the derived class.

	Attributes:
		mesh0 (:class:`dolfin.cpp.mesh.Mesh`): A mesh of the initial condition
		params (:obj:`dict`): A dictionary of parameters
	"""

	def __init__(self, mesh0 = UnitCubeMesh(8,8,8), params={}):
		parameters['form_compiler']['representation'] = 'uflacs'
		parameters['form_compiler']['optimize'] = True
		parameters['form_compiler']['quadrature_degree'] = 4

		self.mesh0 = Mesh(mesh0)
		self.mesh = Mesh(mesh0)

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
			dt (:obj:`float`): A time-step, passed by the :class:`morphosolver.core.simulator.Simulator` object.
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

		ALE.move(self.mesh, self.u)

		self.v.vector()[:] = self.u.vector()[:] / dt
		self.n = FacetNormal(self.mesh)

	def plotMesh(self):
		"""
		Virtual iupdate function. Is run at every time-step of the simulation before solving the growth mechanics\
		by the call to the :func:`~solve` function.
		"""
		plt.figure()
		plot(self.mesh)

# class GrowingTissueBase(TissueBase):
# 	def __init__(self, *args, **kwargs):
# 		super().__init__(*args,**kwargs)

# 		self.a = Constant(0)
# 		self.g = Constant(0)
# 		self.gamma = Constant([0,0,1])

# 	def updateG(self, dt):
# 		g = project(self.g, self.C)
# 		a = project(self.a, self.C)
		
# 		# # print(" ".join( [ str(x) for x in g.vector()[:]] ) )	

# 		gamma = project(self.gamma / sqrt( inner(self.gamma, self.gamma) ), self.V)

# 		self.gFile << (g,self.t)

# 		g0 = g * ( 1 -		a )
# 		g1 = g * ( 1 + 2 *	a )
		
# 		gf = project( Identity(3) * g0 + (g1-g0) * outer( gamma, gamma ), self.W )


# 		self.G = project( (Function(self.W0, gf.vector()) * dt + Identity(3)) * self.G, self.W0 )


# class Tissue(object):
# 	def __init__(self, base):
# 		self.base = base

# 	def __getattr__(self, key):
# 		return getattr(self.base, key)

# 	def __setattr__(self, key, value):
# 		if key == "base":
# 			super().__setattr__(key, value)
# 		else:
# 			setattr(self.base, key, value)

# class OsmoticTissueBase(GrowingTissueBase):
# 	def __init__(self, *args, **kwargs):
# 		super().__init__(*args,**kwargs)

# 		self.c = SmartFunction(self.outputPath, "c", self.C0, functionSpace = self.C)
# 		self.ph = SmartFunction(self.outputPath, "ph", self.C0, functionSpace = self.C) 
# 		self.P = SmartFunction(self.outputPath, "P", self.C0, functionSpace = self.C)
		
# 		self.w = SmartFunction(self.outputPath, "water", self.V0, functionSpace = self.V)

# 		self.g = SmartFunction(self.outputPath, "g", self.C0, functionSpace = self.C)

# 		self.y = SmartFunction(self.outputPath, "y", self.C0, functionSpace = self.C)
		
# 		for param in ["kw", "phi", "rho", "mu", "Dc", "r"]:
# 			if param in self.params:
# 				self.params[param] = SmartFunction(self.outputPath, param, self.C0, function = project(self.params[param], self.C) )
# 			else:
# 				self.params[param] = SmartFunction(self.outputPath, param, self.C0, functionSpace = self.C)

# 		self.params["kw"].project( Constant(1.0) )
# 		self.params["phi"].project( Constant(1.0) )

# 		self.waterPermeable = MeshFunction("size_t", self.mesh, 2)
# 		self.waterPermeable.set_all(0)

# 	def S(self):
# 		C1 = self.params['C1']
# 		return C1*self.A() - C1 * inv(self.A()) - self.ph.initial * Identity(3)

# 	def updateP(self, dt):
# 		u = TrialFunction(self.C)
# 		v = TestFunction(self.C)

# 		f = -self.g

# 		bc = DirichletBC(self.C, Constant(0.0), self.waterPermeable, 1)

# 		a = inner(grad(u), grad(v))*dx
# 		L = f*v*dx

# 		solve(a == L, self.P, bc)

# 		self.ph.project( self.P + self.c )

# 	def updateC(self, dt):
# 		c0 = Function(self.C, Vector(self.c.vector()) )

# 		c = TrialFunction(self.C)
# 		d = TestFunction(self.C)

# 		g = self.g
# 		w = self.w
# 		v = self.v

# 		kw = self.params["kw"]
# 		phi = self.params["phi"]
# 		A = kw*phi / (kw+phi)

# 		D = self.params["Dc"]
# 		rho = self.params["rho"]
# 		mu = self.params["mu"]

# 		r = self.params["r"]

# 		derivativeTerm = (c-c0)*d*dt*dx
# 		advectionTerm = inner( v + mu*w, grad(c0))*d*dt*dx
# 		diffusionTerm = D*inner( grad(c), grad(d) )*dt*dx
# 		growthTerm = div(v)*( rho + c )*d*dt*dx
# 		sourceTerm = -r*d*dt*dx

# 		F = derivativeTerm
# 		# F += advectionTerm
# 		F += diffusionTerm
# 		F += growthTerm
# 		F += sourceTerm

# 		a, L = lhs(F), rhs(F)

# 		solve(a==L, self.c)

# 		self.c.vector()[self.c.vector()<0] = 0

# 	def updateX(self, dt):
# 		self.u = Function(self.V0)
# 		w = TestFunction(self.V0)
# 		du = TrialFunction(self.V0)

# 		x = SpatialCoordinate(self.mesh0)

# 		L = inner( self.S(), self.eps(w) )*dx(degree=4)\
# 		- inner( self.b, w )*dx(degree=4)\
# 		- inner( self.h, w )*ds(degree=4)\
# 		+ inner( 1e-6*self.u, w )*ds(degree=4)\
# 		- inner( Min(x[2]+self.ut[2]+self.u[2], 0) * Constant((0,0,-1.0)), w )*ds(degree=4)

# 		a = derivative(L, self.u, du)

# 		problem = NonlinearVariationalProblem(L, self.u, bcs=[], J=a)
# 		solver = NonlinearVariationalSolver(problem)

# 		solver.solve()

# 		self.ut.vector()[:] = self.ut.vector()[:] + self.u.vector()[:]

# 		self.updateState(dt)

# 	def solve(self, dt):
# 		print("\tUpdating c...")
# 		self.updateC(dt)
		
# 		print("\tUpdating P and w...")
# 		self.updateP(dt)
# 		self.w.project( -grad(self.P) )

# 		print("\tCalculating G...")
# 		kw = self.params["kw"]
# 		phi = self.params["phi"]

# 		A = kw*phi / (kw + phi)

# 		self.g.project( conditional( ge(self.c, 0), -0.33*A*Min( div( grad(self.c) ), 0 ), 0 ) )
# 		self.updateG(dt)

# 		print("\tFinding equilibrium...")
# 		self.updateX(dt)

# 		self.t += dt

