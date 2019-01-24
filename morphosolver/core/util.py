import os
import sys

from dolfin import *

class SmartFunction(Function):
	""" Smart function is a derivative of the :class:`dolfin.cpp.function` with an extended functionality.\
	It allows to specify two function spaces (:class:`dolfin.cpp.function.FunctionSpace`) for both an\
	initial and current conditions. The state of the Smart functions are automatically saved at each simulation step\
	to the :attr:`morphosolver.core.simulator.Simulator.outputDir`
	
	Attributes:
		name (:obj:`str`): Name of the function. The name is used for a filename, where the function state will be saved.
		functionSpace0 (:class:`dolfin.cpp.function.FunctionSpace`): Function space of the initial condition.
		functionSpace (:class:`dolfin.cpp.function.FunctionSpace`): Function space of the current condition. Either supplied at\
		the initialisation or is derived from the supplied :attr:`~function`.
	"""
	def __init__(self, name, functionSpace0 = None, function = None, functionSpace = None):
		
		self.userName = name
		self.file = File(os.path.join( sys.path[0], "{0}.pvd".format(self.userName) ))

		if function:
			super().__init__(function.function_space(), function.vector())
		elif functionSpace:
			super().__init__(functionSpace)
		else:
			raise Exception("Either function or functionSpace must be supplied.")
		
		self.functionSpace0 = None
		if functionSpace0:
			self.functionSpace0 = functionSpace0
			self.initial = Function(functionSpace0, self.vector())

		self.functionSpace = self.function_space()

	def setOutputPath(self, outputPath):
		"""Set outut directory for saving the function state. Usually is done automatically by the simulator.

		Args:
			outputPath(:obj:`str`): A path to the output directory.
		"""
		self.file = File( os.path.join(outputPath, "{0}.pvd".format(self.userName) ))

	def save(self, t):
		""" Save function state.

		Args:
			t(:obj:`float`): Time associated with the current state.
		"""

		self.file << (self, t)

	def project(self, obj, initial = False):
		"""Project given function to the Smart Function.

		Args:
			obj(:class:`dolfin.cpp.function.GenericFunction`): Function to project.
			initial(:obj:`bool`): Project to initial function space. Defaults to False.
		"""
		if initial and self.functionSpace0:
			self.vector()[:] = project(obj, self.functionSpace0).vector()[:]
		else:
			self.vector()[:] = project(obj, self.functionSpace).vector()[:]

class GenericDomain(SubDomain):
	""" A helper class to define a subdomain with a given function.

	Arguments:
		func(:obj:`lambda(x, onBoundary)`): A lambda function, which returns true if the point x is inside the defined subdomain.
	"""
	def __init__(self, func):
		super().__init__()
		self.func = func

	def inside(self, x, onBoundary):
		return self.func(x, onBoundary)

# def diffuse(function, D, dt, v, r=Constant(0), advection=False, dilution=True):
# 	functionSpace = function.function_space()

# 	c0 = Function(functionSpace, Vector(function.vector()) )

# 	c = TrialFunction(functionSpace)
# 	d = TestFunction(functionSpace)

# 	diffusionTerm	=  dt*D*inner( grad(c), grad(d) )*dx
# 	derivativeTerm  =     d*(c-c0)*dx
# 	reactionTerm	= -dt*d*r*dx
# 	advectionTerm	=  dt*d*inner(v, grad(c))*dx
# 	dilutionTerm	=  dt*d*c*div(v)*dx

# 	F = derivativeTerm + diffusionTerm + reactionTerm
	
# 	if advection:
# 		F = F + advectionTerm

# 	if dilution:
# 		F = F + dilutionTerm

# 	a, L = lhs(F), rhs(F)

# 	c = Function(functionSpace)

# 	solve(a==L, c)

# 	c.vector()[c.vector()<0] = 0

# 	function.vector()[:] = c.vector()[:]

# def diffuseRunge(function, D, dt, v, r=0, advection=False, dilution=True):
# 	functionSpace = function.function_space()

# 	c0 = Function(functionSpace, Vector(function.vector()) )

# 	c = TrialFunction(functionSpace)
# 	d = TestFunction(functionSpace)

# 	diffusionTerm	=  dt*D*inner( grad(c), grad(d) )*dx
# 	derivativeTerm  =     d*(c-c0)*dx
# 	reactionTerm	= -dt*d*r*dx
# 	advectionTerm	=  dt*d*inner(v, grad(c))*dx
# 	dilutionTerm	=  dt*d*c*div(v)*dx

# 	diffusionTerm0  =  dt*D*inner( grad(c0), grad(d) )*dx
# 	advectionTerm0  =  dt*d*inner(v, grad(c0))*dx
# 	dilutionTerm0	=  dt*d*c0*div(v)*dx

# 	F = derivativeTerm + reactionTerm + 0.5*( diffusionTerm + diffusionTerm0 )

# 	if advection:
# 		F = F + advectionTerm + advectionTerm0

# 	if dilution:
# 		F = F + dilutionTerm + dilutionTerm0

# 	a, L = lhs(F), rhs(F)

# 	c = Function(functionSpace)

# 	solve(a==L, c)

# 	c.vector()[c.vector()<0] = 0

# 	function.vector()[:] = c.vector()[:]
