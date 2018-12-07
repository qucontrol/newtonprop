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

The Newton propagator evaluates an expansion of the time evolution operator
:math:`e^{-i \Op{H} dt}` or :math:`e^{\Liouville dt}` in Newton polynomials,
using an implicitly restarted Arnoldi scheme. More generally, it can evaluate
the application of any operator-valued function to a state.

Development of the ``newtonprop`` package happens on `Github`_.
You can read the full documentation at `ReadTheDocs`_.

.. Warning::

    This is a reference implementation only. It aims to be easy to understand,
    so that that it can guide the implementation of the algorithm in a compiled
    language, and to provide a baseline against which to test such an
    implementation. Being written in pure Python, it runs several orders of
    magnitude slower than an implementation in a compiled language. Thus, it
    cannot compete e.g. with the `ODE solvers provided by SciPy`_ (which run at
    C speed), even though the Newton propagator is usually expected to be
    superior in runtime, memory usage, and precision.


.. _ReadTheDocs: https://newtonprop.readthedocs.io/en/latest/
.. _ODE solvers provided by SciPy: https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.integrate.ode.html


Prerequisites
-------------

The `newtonprop` only depends on NumPy_. You may consider the use of QuTiP_ for
efficient data structures for quantum states, operators, and super-operators
(see the Example). As an optional dependency, having Numba_ installed when you
import `newtonprop` can considerably speed up the propagation, assuming the
application of the Hamiltonian/Liouvillian is implemented efficiently. Even
then, though, the implementation will not approach C speed.

.. _NumPy: http://www.numpy.org
.. _Numba: http://numba.pydata.org
.. _QuTiP: http://qutip.org


Installation
------------

To install the latest released version of ``newtonprop``, run this command in your terminal:

.. code-block:: console

    $ pip install newtonprop

This is the preferred method to install ``newtonprop``, as it will always install the most recent stable release.

If you don't have `pip`_ installed, the `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


To install the latest development version of ``newtonprop`` from `Github`_:

.. code-block:: console

    $ pip install git+https://github.com/qucontrol/newtonprop.git@master#egg=newtonprop



.. _Github: https://github.com/qucontrol/newtonprop

Usage
-----

The ``newtonprop`` package exposes its functionality through a single function,
accessible either as ``newtonprop.newton`` or ``newtonprop.propagator.step``::

    >>> from newtonprop import newton
