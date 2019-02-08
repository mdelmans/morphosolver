from morphosolver import *
from dolfin import *

import argparse

class AnisotropicGrowthTissue(TissueBase):
	def init(self):
		self.chemical = Chemical(self, 'morphogen')
		self.chemical.setDiffusion(self.params['D'])
		self.chemical.setReaction( Expression("exp(x[2])", degree = 3) - self.params['d'] * self.chemical.c0 )

		self.growthFunction = AnisoGrowthFunction(self)

	def update(self, dt):
		self.chemical.update(dt)

		self.growthFunction.projectAGG(a = Constant(self.params['a']), g = 0.1*self.chemical, gamma = Constant(self.params['g']) )
		self.growthFunction.updateG(dt)

if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description=
	'''
	Morphosolver demo: Anisotropic growth\n\
	-------------------------------------\n\n\
	Simulates reaction diffusion of a chemical [c] with reaction term r = e^(z) - 0.1*c on a UnitCubeMesh\
	''', formatter_class = argparse.RawTextHelpFormatter)

	parser.add_argument('-a',		metavar = 'anisotropy',	type = float,	default = 0.8,	help = 'Anisotropy of a tissue [0,1]. Default a = 0.8' )
	parser.add_argument('-g', metavar=('gx', 'gy', 'gz'), nargs=3, default = [0,0,1], help='Vector of anispotropy. Default g = [0.0, 0.0, 1.0]')

	parser.add_argument('-D',		metavar = 'diffusion',	type = float,	default = 0.1,	help = 'Diffusion rate. Default D = 0.1' )
	parser.add_argument('-d',		metavar = 'deg',		type = float,	default = 0.1,	help = 'Degradaion rate. Default d = 0.1' )

	parser.add_argument('-n',		metavar = 'steps',		type = int,		default = 20, 	help = 'Number of timesteps. Default steps = 20')
	parser.add_argument('-dt',								type = float,	default = 0.1, 	help = 'Timestep. Default dt = 0.1')

	params =  vars(parser.parse_args())

	tissue = AnisotropicGrowthTissue(params = params)
	
	simulator = Simulator(tissue)
	simulator.run(params['dt'], params['n'])