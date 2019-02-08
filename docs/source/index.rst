.. morphosolver documentation master file, created by
   sphinx-quickstart on Wed Jan 23 17:01:01 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Morphosolver
============

A FEniCS-based library for solving reaction-diffusion equations on growing domains. Developed as part of a PhD thesis by Mihails Delmans.

Installation
------------

Morphosolver is build on top of the `FEniCS
<https://fenicsproject.org/>`_ Docker container.

In order to use morphosolver, first install `Docker
<https://www.docker.com/>`_, and then run the following:

.. code-block:: console
   
   $ docker run -it -v /path/to/home:/morphosolverhome mdelmans/morphosolver:2018.1

Examples
--------

See examples in `Models section <models/models.html>`__.

.. toctree::
   :maxdepth: 3
   :caption: Content

   installation
   modules/morphosolver
   models/models



