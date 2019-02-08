from ...core.util import SmartFunction
from .base import GrowthFunctionBase

from dolfin import Constant, Identity
from dolfin import outer, inner, sqrt

class AnisoGrowthFunction(GrowthFunctionBase):
	""" Anisotropic growth function.
		:math:`\mathcal{G} = g_0\mathbf{I} + (g_1-g_0)\gamma \otimes \gamma`
		
		where :math:`g_0 = g(1-a)`, :math:`g_1 = g(1+2a)`
		
		:math:`3g` is relative volumetric expansion
		
		:math:`a` is anisotropy constant in range [0,1]
		
		:math:`\gamma` is an anisotropy vector

		Args:
			tissue (:class:`~morphosolver.core.tissuebase.TissueBase`): A tissue to operate on
	"""
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
		""" Calculate and project growth function from given :math:`a`, :math:`g` and :math:`\gamma` values
		"""
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