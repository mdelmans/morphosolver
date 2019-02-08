AnisotropicGrowth
=================

Problem
-------

Consider a unit cube tissue, where diffusible growth-inducing morphogen is produced with a rate proportional to the z-axis:

:math:`\frac{dc}{dt} = D\nabla^2c + e^z - dc`

where :math:`c` is concentration of a morphogen, :math:`D` is a diffusion constant and :math:`d` is a degradation rate.

Initialise
----------

To implement this problem using `morphosolver`, we should create a new tissue class (:class:`~morphosolver.core.tissuebase.TissueBase`) add a new chemical (:class:`~morphosolver.core.chemical.Chemical`) and set it's diffusion and reaction rates using :func:`~morphosolver.core.chemical.Chemical.setDiffusion` and :func:`~morphosolver.core.chemical.Chemical.setReaction` methods during the initialisation stage:

.. literalinclude:: ../../../models/anisotropicgrowth/model.py
	:lines: 6-10

We then add an anisotropic growth function(:class:`~morphosolver.addons.growthfunctions.AnisoGrowthFunction`) object to the tissue:

.. literalinclude:: ../../../models/anisotropicgrowth/model.py
	:lines: 12

Update
------

Then during the update step, we should update the chemical state using the :func:`~morphosolver.core.chemical.Chemical.update` method:

.. literalinclude:: ../../../models/anisotropicgrowth/model.py
	:lines: 14-15

apply growth rate proportional to the chemical concentration using :func:`~morphosolver.addons.growthfunctions.AnisoGrowthFunction.projectAGG` method and update the growth tensor :math:`\mathbf{G}` of the tissue using the :func:`~morphosolver.addons.growthfunctions.GrowthFunctionBase.updateG` method of the growth function:

.. literalinclude:: ../../../models/anisotropicgrowth/model.py
	:lines: 17-18

Run
---

To simulate the problem, we create a new instance of the :class:`AnisotropicGrowthTissue`:

.. literalinclude:: ../../../models/anisotropicgrowth/model.py
	:lines: 40

create and run the simulation:

.. literalinclude:: ../../../models/anisotropicgrowth/model.py
	:lines: 42-43

Complete listing
----------------

.. literalinclude:: ../../../models/anisotropicgrowth/model.py
	:linenos: