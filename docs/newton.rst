The Newton Propagator
=====================

*The following overview has been adapted from Ref* :cite:`GoerzPhd2015`.

Quantum mechanical equations of motion
--------------------------------------

The behavior of a quantum system is usually described either by the Schrödinger equation

.. math::

    \frac{\partial}{\partial t} \ket{\Psi(t)} = -\frac{\ii}{\hbar} \Op{H} \ket{\Psi(t)},

or, for open quantum systems, by the Liouville-von-Neumann equation

.. math::

    \frac{\partial}{\partial t} \Op{\rho}(t) = \Liouville \Op{\rho}(t).


ODE solvers versus series expansion
-----------------------------------

There are two possible approaches to obtaining a solution. The first is
to simply take the equation of motion and apply one of the generic
numerical methods for solving ordinary differential equations (ODEs),
like one of the Runge-Kutta (RK)
methods :cite:`LambertODEBook,NumRecipesFortran`.

The alternative approach is to solve the equation of motion
analytically, and then to evaluate that solution numerically. In this
way, results of arbitrary precision can be obtained, with the obvious
caveat that the propagation scheme will be specific to a particular
equation of motion and its solution.

The Schrödinger equation for a time-independent Hamiltonian has the
solution

.. math::

   \Ket{\Psi(t)}= \Op{U}(t,0) \Ket{\Psi(0)}
                  = \ee^{-\frac{\ii}{\hbar} \Op{H}\, t} \Ket{\Psi(0)}\,.

Similarly, the Liouville-von-Neumann equation has the solution

.. math::

   \Op{\rho}(t) = \DynMap(t,0) \
                  = \ee^{\Liouville\, t} \Op{\rho}(0)\,.


For time-dependent Hamiltonians, we approximate :math:`\Op{H}(t)` as
piecewise constant on a time grid with time step :math:`dt`. Then, the
total time evolution operator is the product of the time evolution
operators at each time step,

.. math::

   \Op{U}[T,0]
     = \prod_{i=1}^{nt-1} \underbrace{\Op{U}(t_i+dt, t_i)}_{\equiv \Op{U}_i}
     = \prod_{i=1}^{nt-1} \exp\bigg[-\frac{\ii}{\hbar}
       \underbrace{\Op{H}\left(t_i + \frac{dt}{2}\right)}_{\equiv \Op{H}_i} dt\bigg]\,.

The time step :math:`dt` must be chosen sufficiently small that this is
a good approximation; in practice, convergence is checked by verifying
that the numerical results remain stable within a desired precision when
:math:`dt` is decreased. The same applies to the Liouvillian.

For a system of non-trivial size,
:math:`\exp[-\frac{\ii}{\hbar} \Op{H} dt]` is evaluated by expanding the
exponential in a polynomial series,

.. math::

   \ee^{-\frac{\ii}{\hbar} \Op{H}\, dt} \Ket{\Psi}
     \approx \sum_{n=0}^{N-1} a_n P_n(\Op{H}) \Ket{\Psi}\,.
     \label{eq:poly_expansion}

where :math:`P_n(\Op{H})` is a polynomial of degree :math:`n` and
:math:`\{a_n\}` are the expansion coefficients. Applying
:math:`P_n(\Op{H})` to :math:`\ket{\Psi}` then simply means repeated
applications of :math:`\Op{H}`. For this reason, an efficient
propagation relies on the proper use of sparsity in storing the
Hamiltonian and spectral methods such as the Fourier grid.

The idea of evaluating the exponential as a polynomial series is already
presupposed by the very definition of the exponential of an operator,

.. math:: \exp[\Op{A}] \equiv \sum_{n=0}^{\infty} \frac{1}{n!} \Op{A}^n\,.

However, this expansion converges particularly slowly and is
numerically unstable :cite:`Tal-EzerJCP84`. Thus, it is not
suitable for time propagation. Instead, a polynomial basis must be
chosen such that the series converges quickly and can be truncated as early as
possible.


Chebychev expansion
-------------------

For a function :math:`f(x)` with :math:`x \in [-1, 1]`, it can be
shown :cite:`GilBook2007` that the fastest converging
polynomial series is the expansion in Chebychev polynomials.
When using the Chebychev expansion for propagation, the requirement that
the argument of :math:`f(x)` must be real translates into :math:`\Op{H}`
being Hermitian. Secondly, to account for the requirement that
:math:`x \in [-1, 1]`, the Hamiltonian must be normalized
as :cite:`KosloffJCP88,TannorBook,NdongJCP09`

.. math:: \Op{H}_{\text{norm}} = 2 \frac{\Op{H} - E_{\min}\,\identity}{\Delta} - \identity\,,

where :math:`\Delta = E_{\max} - E_{\min}` is the spectral radius and
:math:`E_{\max}` and :math:`E_{\min}` are the smallest and largest
eigenvalue.

Newton expansion
----------------

For a general function :math:`f(z)` with :math:`z \in \Complex`, expansion into
Chebychev polynomial is not possible. For a function defined on an operator,
a complex :math:`z` corresponds to an operator with complex eigenvalues.
This is the case for a non-Hermitian Hamiltonian, and, more importantly, for a
dissipative Liouvillian in the Liouville-von-Neumann equation. In this case,
Newton polynomials must be used. The expansion in Newton polynomials
:math:`R_n(z)` reads

.. math::

   f(z) \approx \sum_{n=0}^{N-1} a_n R_n(z)\,, \quad
     R_n(z) = \prod_{j=0}^{n-1} \left( z-z_j \right)\,,

for a set of sampling points :math:`\{z_j\}` at which the interpolation
is exact. The coefficients are defined as the divided differences :cite:`AshkenaziJCP95`,

.. math::
   :label: divided_differences

   \begin{aligned}
     a_0 &= f(z_0)\,, \\
     a_1 &= f(z_1) - f(z_0)\,, \\
     a_n &= \frac{f(z_n) - \sum_{j=0}^{n-1} a_j
                  \prod_{k=0}^{j-1}\left(z_n - z_k\right)}
                 {\prod_{j=0}^{n-1} \left(z_n - z_j\right)}\,.\end{aligned}

For solving the Liouville-von Neumann equation,
:math:`f(z)=\ee^{z \,dt}`, where the argument :math:`z` is
:math:`\Liouville`. Thus, the propagation is written as

.. math::

   \Op{\rho} = \ee^{\Liouville \,dt} \, \Op{\rho}_0
               \approx
                 \underbrace{%
                   \sum_{n=0}^{N-1} a_n
                   \left( \Liouville - z_n\identity \right) }_{%
                             \equiv p_{N-1}(\Liouville)}
                   \Op{\rho}_0\,,

where the polynomial is evaluated through repeated application the Liouvillian.

The central issue for obtaining a fast-converging series is a proper
choice of the sampling points :math:`\{z_j\}`. The fastest convergence
results from using the complex eigenvalues of
:math:`\Liouville` :cite:`KosloffARPC94`. However, the exact
eigenvalues of the Liouvillian are not readily available. More
generally, arbitrary points from the spectral domain of
:math:`\Liouville` can be used as sampling points.

A widely used method is to estimate the spectral domain and to encircle
it with a rectangle or ellipse :cite:`BermanJPA92,AshkenaziJCP95,HuisingaJCP99`.
Then, a large number of expansion coefficients are
calculated from sampling points on that boundary. The same coefficients
are used for the propagation of all Liouvillians on the time grid under
the assumption that they all fit into the same encirclement. The series
is truncated as soon as convergence is reached. This is similar to the
method employed for the Chebychev propagator, where a set of
coefficients is calculated once and then used for the propagation of any
Hamiltonian that is within the same spectral range.

A middle path between the exact eigenvalues of :math:`\Liouville` and
the crude encirclement of the spectral domain is the use of the Krylov
method to obtain approximate eigenvalues.
The `Arnoldi algorithm`_ for
:math:`\hat{A} = \Liouville` and using :math:`\vec{v} = \Op{\rho}` as a starting vector
yields a set of approximate eigenvalues of :math:`\Liouville`, as well
as a Hessenberg matrix :math:`\hat{H}` that is the projection of
:math:`\Liouville` into the Krylov subspace, and the set of Arnoldi
vectors that span that subspace. Instead of using :math:`\Liouville` as
the argument of the polynomial :math:`p_{N-1}`, the Hessenberg matrix
may be used. If :math:`\hat{V}_{N}` is the transformation matrix between
the full Liouville space and the reduced Krylov space, consisting of the
Arnoldi vectors as columns, the propagation is evaluated using

.. math::

   \Liouville \Op{\rho}
     \approx
     \hat{V}_{N}\, p_{m-1}(\hat{H}) \, \hat{V}_{N}^{\dagger} \Op{\rho}_0\,.

Assuming :math:`N` is much smaller than the full dimension of the
Liouville space, most of the numerical effort is in the Arnoldi
algorithm, in constructing the Krylov space.

.. _Arnoldi algorithm: https://en.m.wikipedia.org/wiki/Arnoldi_iteration


Newton propagation with an implicitly restarted Arnoldi method
--------------------------------------------------------------

However, even for moderate values of :math:`N` (typically on the order
of 100), the Arnoldi algorithm can require prohibitive amounts of
memory. This is because a full set of :math:`N` Arnoldi vectors, each of
the dimension of the Liouville space, needs to be stored. To counter this
problem, an iterative scheme has been
developed :cite:`Tal-EzerSJSC2007`. Instead of performing
the Arnoldi algorithm to a high order :math:`N`, until convergence is
reached in the propagation, we stop at some small order :math:`m<10`.
This gives a first approximation to the propagated density matrix,

.. math::
   :label: newton1stIter

   \Op{\rho}^{(1)}
     = p_{m-1}^{(0)}(\Liouville) \Op{\rho}_0
     = \sum_{n=0}^{m-1} a_n R_n(\Liouville) \Op{\rho}_0\,.
     \label{eq:newton1stIter}

The idea is now to iteratively add remaining terms to the Newton series
in chunks of size :math:`m`, retaining all coefficients and sampling
points, but restarting the Arnoldi procedure in every iteration.

Adding the next :math:`m` terms to Eq. :eq:`newton1stIter` yields

.. math::

   \begin{split}
     \Op{\rho}^{(2)}
    &= \Op{\rho}^{(1)} + \sum_{n=m}^{2m-1} a_n R_n(\Liouville) \Op{\rho}_{0} \\
    &= \Op{\rho}^{(1)}
       + \underbrace{\left(\sum_{n=0}^{m-1}  a_{m+n} R_n^{(1)}(\Liouville)\right)}_{%
                                  \equiv p_{m-1}^{(1)} }
         \underbrace{\left(R_{m}^{(0)}(\Liouville) \Op{\rho}_{0} \right)}_{%
                              \equiv \Op{\sigma}^{(1)} }\,,
   \end{split}

with

.. math::

   R_n^{(0)}(\Liouville) = \prod_{j=0}^{n-1}(\Liouville - z_{j}\identity), \qquad
     R_n^{(1)}(\Liouville) = \prod_{j=0}^{n-1}(\Liouville - z_{n+j}\identity)\,.

That is, the terms in :math:`R_n(\Liouville)` already known from the
calculation of :math:`\Op{\rho}^{(1)}` have been pulled out, and yield a
new "starting vector" :math:`\Op{\sigma}^{(1)}`, which is the argument
to a Newton series of only :math:`m` new terms. The new sampling points
on which the :math:`R_n^{(1)}` are evaluated are obtained by applying
the Arnoldi procedure to :math:`\Op{\sigma}^{(1)}`. The Newton
coefficients continue recursively from the previous restart. The third
iteration yields

.. math::

   \Op{\rho}^{(3)}
     = \Op{\rho}^{(2)}
       + \underbrace{\left(\sum_{n=0}^{m-1}  a_{2m+n} R_n^{(2)}(\Liouville)\right)}_{%
                                  \equiv p_{m-1}^{(2)} }
         \underbrace{\left(R_{m}^{(1)}(\Liouville) \Op{\sigma}_{1} \right)}_{%
                              \equiv \Op{\sigma}^{(2)} }\,.

The Newton propagator continues, adding the :math:`m` terms evaluating

.. math::

   p_{m-1}^{(s)}(\Liouville) \Op{\sigma}^{(s)}
     = \sum_{n=0}^{m-1} a_{sm + n}
       \prod_{k=0}^{n-1} \left(\Liouville - z_{sm+k} \identity \right)
       \Op{\sigma}^{(s)}

with

.. math:: \Op{\sigma}^{(s)} = p_{m-1}^{(s-1)} \Op{\sigma}^{(s-1)}

at every restart iteration.

In the implementation of the algorithm, there are two details that need
to be taken into account for numerical stability. First, the denominator
of the divided differences in Eq. :eq:`divided_differences` may become extremely small if
consecutive sampling points are close to each other. This can be
addressed by reordering the points such that the denominator in the
divided differences is maximized. This process is called Leja
ordering :cite:`ReichelBIT1990`. The reverse problem that
the sampling points are too far apart, causing an underflow in the
calculation of coefficients can be avoided by normalizing the
Liouvillian as

.. math:: \tilde{\Liouville} = \frac{1}{\rho} \left( \Liouville - c \right)\,,

where :math:`c` is an estimate for the center of the spectrum of
:math:`\Liouville`, and the eigenvalues are roughly contained in a
radius :math:`\rho` around :math:`c`. These values can be estimated from
the sampling points obtained in the first iteration of the Newton
propagator. The normalization of the Liouvillian is in some sense
similar to the normalization of the Hamiltonian in the Chebychev
propagator, but it is crucial there since the Chebychev polynomials are
only defined in the domain :math:`[-1, 1]`. For the Newton propagator,
the normalization is only for numerical stability.

References
----------

.. bibliography:: refs.bib
   :cited:
   :style: unsrt
