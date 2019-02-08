from dolfin import Constant, Identity

from ...core.util import SmartFunction

class GrowthFunctionBase(SmartFunction):
	""" Base of the growth function (:math:`\mathcal{G}`). Operates on the groth tensor `G` of a tissue, such that:
		:math:`\mathbf{G}(t+\delta t) = \mathbf{G}(t)(\mathcal{G}(t) dt + \mathbf{I})`
		
		Args:
			tissue (:class:`~morphosolver.core.tissuebase.TissueBase`): A tissue to operate on
	"""
	def __init__(self, tissue):
		self.tissue = tissue
		super().__init__("Gf", functionSpace = self.tissue.W, functionSpace0 = self.tissue.W0)

	def updateG(self, dt):
		""" Update the G tensor of a tissue.

		Args:
			dt (:obj:`float`): Timestep
		"""
		self.tissue.G.project( (self.initial * dt + Identity(3)) * self.tissue.G )