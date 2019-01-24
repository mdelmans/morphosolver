from dolfin import *
from .util import SmartFunction


class Chemical(SmartFunction):
	""" Chemical defines a concentration field which can be subject to diffusion and user-specified reaction.
	
	Attributes:
		tissue (:class:`morphosolver.core.tissuebase.TissueBase`): the host
		name (:obj:`str`): Name of the chemical
	"""
	def __init__(self, tissue, name):
		super().__init__( "ch_{0}".format(name), functionSpace = tissue.C )

		self.c0 = Function(self.functionSpace)
		self.setDiffusion(0)
		self.setReaction(0)
		self.setExtra(0)

	def setDiffusion(self, D):
		"""Set diffusion

		Args:
			D (:obj:`float`): Diffusion rate
		"""
		self.D = D

	def setReaction(self, R):
		"""Set reaction term
		
		Args:
			R (:class:`dolfin.cpp.function.GenericFunction`): Reaction term
		"""
		self.R = R

	def setExtra(self, E):
		self.E = E

	def update(self, dt):
		""" Update function by calculating reaction-diffusion step

		Args:
			dt (:obj:`float`): Time-step
		"""
		self.c0.vector()[:] = self.vector()[:]

		c  = TrialFunction(self.functionSpace)
		e  = TestFunction(self.functionSpace) 

		derivativeTerm = e*(c-self.c0)*dx

		F = derivativeTerm
		
		try:
			diffusionTerm  = dt*self.D*inner( grad(c), grad(e) )*dx
			F = F + diffusionTerm
		except Exception as ex:
			print("Diffusion term: ", ex)
		
		try:
			reactionTerm   = -dt*self.R*e*dx
			F = F + reactionTerm
		except Exception as ex:
			print("Reaction term: ", ex)

		if self.E:
			F = F + E(c,e)

		a, L = lhs(F), rhs(F)

		c = Function(self.functionSpace)

		solve(a==L, c)

		c.vector()[c.vector() < 0] = 0
		
		self.vector()[:] = c.vector()[:]