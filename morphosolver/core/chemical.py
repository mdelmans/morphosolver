from dolfin import *
from .util import SmartFunction


class Chemical(SmartFunction):
	def __init__(self, tissue, name):
		super().__init__( "ch_{0}".format(name), functionSpace = tissue.C )

		self.c0 = Function(self.functionSpace)
		self.setDiffusion(0)
		self.setReaction(0)
		self.setExtra(0)

	def setDiffusion(self, D):
		self.D = D

	def setReaction(self, R):
		self.R = R

	def setExtra(self, E):
		self.E = E

	def update(self, dt):
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

# class Chemical: 
# 	def __init__(self, tissue, name):
# 		self.name	= name
# 		self.tissue	= tissue
# 		self.c		= SmartFunction(self.tissue.outputPath, "ci_{0}".format(name), functionSpace = self.tissue.C)
# 		self.c0		= Function(self.tissue.C)

# 		self.setDiffusion(0)
# 		self.setReaction(0)

# 	def setDiffusion(self, D):
# 		self.D = D

# 	def setReaction(self, R):
# 		self.R = R

# 	def update(self, dt):
# 		self.c0.vector()[:] = self.c.vector()[:]

# 		c  = TrialFunction(self.tissue.C)
# 		e  = TestFunction(self.tissue.C) 

# 		derivativeTerm = e*(c-self.c0)*dx

# 		F = derivativeTerm
		
# 		try:
# 			diffusionTerm  = dt*self.D*inner( grad(c), grad(e) )*dx
# 			F = F + diffusionTerm
# 		except Exception as ex:
# 			print("Diffusion term: ", ex)
		
# 		try:
# 			reactionTerm   = -dt*self.R*e*dx
# 			F = F + reactionTerm
# 		except Exception as ex:
# 			print("Reaction term: ", ex)

# 		a, L = lhs(F), rhs(F)

# 		c = Function(self.tissue.C)

# 		solve(a==L, c)

# 		c.vector()[c.vector() < 0] = 0
		
# 		self.c.vector()[:] = c.vector()[:]

# 		self.c.save(self.tissue.t)

# 	def project(self, c):
# 		self.c.project(c)