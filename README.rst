=============================
The newtonprop Python package
=============================

.. image:: https://img.shields.io/badge/github-qucontrol/newtonprop-blue.svg
   :alt: Source code on Github
   :target: https://github.com/qucontrol/newtonprop
.. image:: https://img.shields.io/pypi/v/newtonprop.svg
   :alt: newtonprop on the Python Package Index
   :target: https://pypi.python.org/pypi/newtonprop
.. image:: https://img.shields.io/travis/qucontrol/newtonprop.svg
   :alt: Travis Continuous Integration
   :target: https://travis-ci.org/qucontrol/newtonprop
.. image:: https://ci.appveyor.com/api/projects/status/vf9vb3k6dqee1oad?svg=true
   :alt: AppVeyor Continuous Integration
   :target: https://ci.appveyor.com/project/goerz/newtonprop
.. image:: https://img.shields.io/coveralls/github/qucontrol/newtonprop/master.svg
   :alt: Coveralls
   :target: https://coveralls.io/github/qucontrol/newtonprop?branch=master
.. image:: https://readthedocs.org/projects/newtonprop/badge/?version=latest
   :alt: Documentation Status
   :target: https://newtonprop.readthedocs.io/en/latest/?badge=latest
.. image:: https://img.shields.io/badge/License-BSD-green.svg
   :alt: BSD License
   :target: https://opensource.org/licenses/BSD-3-Clause

Pure Python reference implementation of the Newton propagator for quantum dynamics.

Development of `newtonprop` happens on `Github`_.
You can read the full documentation at `ReadTheDocs`_.


.. _ReadTheDocs: https://newtonprop.readthedocs.io/en/latest/


Installation
------------
To install the latest released version of ``newtonprop``, run this command in your terminal:

.. code-block:: console

    $ pip install newtonprop

This is the preferred method to install ``newtonprop``, as it will always install the most recent stable release.

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


To install the latest development version of ``newtonprop`` from `Github`_.

.. code-block:: console

    $ pip install git+https://github.com/qucontrol/newtonprop.git@master#egg=newtonprop



.. _Github: https://github.com/qucontrol/newtonprop

Usage
-----

The ``newtonprop`` package exposes its functionality through a single function,
accessible either as ``newtonprop.newton`` or ``newtonprop.propagator.step``::

    >>> from newtonprop import newton
