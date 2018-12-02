"""Newton Propagator"""
import logging
import numpy as np

__all__ = ['step']


def _arnoldi(A, dt, v0, m_max, inner=np.vdot, norm=np.linalg.norm):
    """Calculate the (extended) Hessenberg matrix of an operator A(t) and
    return it together with the  m+1 (m <= m_max) Arnoldi vectors
    (orthonormlized Krylov vectors).  Also return the combined Ritz values (all
    eigenvalues of the Hessenberg matrices of size 1..m )

    Args
        A (callable): Function encoding the operator A. Calling `A(v0)`
            must return the result of applying A to v0.
        dt (float): Time step
        v0: Initial state
        m_max (int): Maximal Krylov dimension
        inner (callable): Function that evaluates an inner product. Must take
            two arguments of the type of `v0` and return a complex number. If
            None, defaults to :func:`numpy.vdot`.
        norm (callable): Function that calculates the norm of a state. Must
            take one argument of the type of `v0` and return a real number. If
            None, defaults to :func:`numpy.linalg.norm`.
    """
    m = m_max

    # Hessenberg matrix (at maximum size)
    Hess = np.matrix(np.zeros(shape=(m_max+1, m_max+1), dtype=np.complex128))

    # Eigenvalues of all Hess
    Ritz = []

    arnoldi_vecs = []

    beta = norm(v0)
    if (abs(beta-1.0) > 1.0e-10):
        print("beta = ", beta)
        raise AssertionError(
            "v0 must have norm 1.0. Mismatch between `inner` and `norm`?")
    v = v0 / beta
    arnoldi_vecs.append(v)
    for j in range(m):
        v = A(v)  # v_{j+1}
        for i, v_i in enumerate(arnoldi_vecs):
            Hess[i, j] = dt * inner(v_i, v)
            v = v - (Hess[i, j]/dt) * v_i
        # At this point, we have finished the (j+1) x (j+1) Hessenberg matrix
        Ritz.extend(np.linalg.eigvals(Hess[:j+1, :j+1]))
        h_next = norm(v)
        Hess[j+1, j] = h_next * dt
        if h_next <= 1.0e-14:  # abort early at convergence
            m = j
            break
        v *= 1 / h_next  # normalize
        arnoldi_vecs.append(v)
    # At this point, arnoldi_vecs contains m+1 elements
    Ritz = np.array(Ritz, dtype=np.complex128)
    return arnoldi_vecs, Hess[:m+1, :m+1], Ritz, m


def _normalize_points(z):
    """Given a set of complex points z, return the normalization radius"""
    r = 0.0
    for z_i in z:
        r_i = abs(z_i)
        if r_i > r:
            r = r_i
    # we need to enlarge the radius a little bit to account for points that
    # will be added in later iterations
    r *= 1.2  # arbitary factor
    assert(r > 0.0), "Radius is zero"
    return r


def _extend_newton_coeffs(func, old_a, new_leja, center, radius):
    """Extend a set of Newton coefficients, by using a set of new_leja points
    which are normalized with the given center and radius
    """
    n_old = len(old_a)
    m = len(new_leja) - n_old
    a = np.zeros(n_old+m, dtype=np.complex128)
    a[:n_old] = old_a
    n0 = n_old

    if n_old == 0:
        a[0] = func(new_leja[0])
        n0 = 1

    for k in range(n0, n_old+m):
        d = 1.0
        pn = 0.0
        for n in range(1, k):  # 1..k-1
            zd = new_leja[k] - new_leja[n-1]
            d *= zd / radius
            pn += a[n] * d
        zd = new_leja[k] - new_leja[k-1]
        d *= zd / radius
        assert(abs(d) > 1.0e-200), "Divided differences too small"
        a[k] = (func(new_leja[k]) - a[0] - pn) / d
    return a


def _extend_leja(old_leja, new_points, n_use):
    """Given a set of normalized (ordered) Leja points, extract `n_use` points
    from the (normalized) new_points, and append them to the set of leja points
    """
    n_old = len(old_leja)
    new_leja = np.zeros(n_old + n_use, dtype=np.complex128)
    new_leja[:n_old] = old_leja[:]
    i_add_start = 0
    if n_old == 0:
        # At the very beginning, start with the point that has largest absolute
        # value
        for i in range(len(new_points)-1):  # 0 .. n_old - 2
            if (abs(new_points[i]) > abs(new_points[-1])):
                temp = new_points[i]
                new_points[i] = new_points[-1]
                new_points[-1] = temp
        new_leja[0] = new_points[-1]
        i_add_start = 1
    # find the best point for new_leja[n_old+n]
    n_added = i_add_start
    ex = 1.0/(n_old + n_use)
    for i_add in range(i_add_start, n_use):
        p_max = 0.0
        i_max = 0
        # the new leja are defined with index  0 .. (n_old-1)+n
        # the new candidates are defined with index 0 .. len(new_points)-1+n
        for i in range(len(new_points)-i_add):  # trial points (candidates)
            p = 1.0
            for j in range(n_old + i_add):  # existing leja points
                p *= np.abs(new_points[i] - new_leja[j])**ex
            # at this point p is the divided difference denominator for the
            # candidate with index i
            if p > p_max:
                p_max = p
                i_max = i
        # XXX if p_max below limit: abort
        new_leja[n_old+i_add] = new_points[i_max]
        n_added += 1
        # remove the used point by moving in the last point
        new_points[i_max] = new_points[len(new_points)-1-i_add]
    return new_leja, n_added


def step(
        A, v, dt, func=None, m_max=10, maxrestart=100, tol=1e-12,
        inner=None, norm=None, zero=None):
    r"""Evaluate $f(A\,dt)\,v$, for example $e^{-i\,A\,dt}\,v$ for the
    propagation of a quantum state.

    Applies the result of an arbitrary abstract operator function $f(A\,dt)$ to
    an abstract state $v$, where $A$ is an abstract operator and $dt$ is a
    scalar factor. The intended application is computing a single time step
    $dt$ of the time evolution of a quantum state, where the form of the time
    evolution operator/dynamical map is $e^{-i\,A\,dt}$. If $v$ is a Hilbert
    space state, $A$ is the system Hamiltonian (constant over the duration
    :math:`dt`). If $v$ is a density matrix, $A = i \Liouville$ with the
    Liouvillian super-operator $\Liouville$.

    The operator $f(A\,dt)$ is never explicitly constructed, only its
    application to $v$ is calculated. This distinguishes the Newton propagator
    from functions like :func:`scipy.linalg.expm` for matrix exponentiation.
    Internally, it uses an iterative restarted Krylov projection of $(A\,dt)$,
    where the Krylov dimension in each iteration is `m_max`. Consequently,
    storage for only `m_max` Arnoldi-vectors of the same size as `v` is
    required.

    Args:
        A (callable): Function encoding the abstract operator $A$. Calling
            ``A(v)`` must return the result of applying $A$ to $v$.
        v: Initial state $v$.
        dt (float): scalar parameter (time step :math:`dt`)
        func (callable): The scalar version of the operator-function $f$ to be
            evaluated, and the result applied to state. Must take one complex
            argument and return a complex value. If None, defaults to
            ``lambda x: numpy.exp(-1j * x)``, that is the function for the
            quantum-mechanical time evolution operator. Note that `func` will
            only ever be called with a scalar argument, not with the argument
            $(A\,dt)$.
        m_max (int): Maximal Krylov dimension.
        maxrestart (int): Maximal number of Newton restarts (iterations of the
            algorithm)
        tol (float): Desired precision of the result.
        inner (callable): Function that evaluates an inner product. Must take
            two arguments of the type of `v0` and return a complex number. If
            None, defaults to :func:`numpy.vdot`.
        norm (callable): Function that calculates the norm of a state. Must
            take one argument of the type of `v0` and return a real number. If
            None, defaults to :func:`numpy.linalg.norm`. The norm must be
            induced by the inner product specified by `inner`.
        zero (callable): Function that takes `v` as input and allocates and
            returns a new zero state of the same type and shape as `v`.
            If None, defaults to
            ``lambda v: numpy.zeros(shape=v.shape, dtype=v.dtype)``

    Returns:
        The result of the operator-function $f(A\,dt)$ applied to $v$.
        The return value will be of the same type as the input `v`.

    .. Note::

        The input `state` may be of any array-like type. It should have (but is
        not strictly required to have) a `shape` attribute. It must support
        being multiplied or divided by a scalar, e.g. ``state / 2.0`` or
        ``2 * state``. Lastly, it must support in-place addition with another
        state, ``state += state2``.  All other mathematical properties of
        `state` can be defined via custom `inner`, `norm`, and `zero`
        functions.

        Mathematically, `state` must be an element a Hilbert space (a "complete
        inner product space"). This includes density matrices, which can be
        interpreted as elements of a Hilbert space provided one chooses an
        appropriate norm and inner product. The norm *must* be the one
        induced by the inner product,

        .. math::

            \Norm{v} \equiv \sqrt{\AbsSq{\Braket{v}{v}}}

        The parameters `inner` and `norm` must fulfill this definition.
        For density matrices, they should be the Hilbert-Schmidt product and
        the Hilbert-Schmidt norm. The "operator norm" has no associated inner
        product. If the operator norm is used, density matrices are not
        elements of a Hilbert space, but of a C* algebra, which would not allow
        for the evaluation of the propagator.
    """

    logger = logging.getLogger('newton')

    if inner is None:
        inner = np.vdot
    if norm is None:
        norm = np.linalg.norm
    if func is None:
        func = lambda x: np.exp(-1j * x)
    if zero is None:
        zero = lambda v: np.zeros(shape=v.shape, dtype=v.dtype)

    try:
        N = np.prod(v.shape)
        assert(m_max <= N), "m_max must be smaller than the system dimension"
    except AttributeError:
        pass

    w = zero(v)                                    # result vector
    Z = np.zeros(0, dtype=np.complex128)               # Leja points
    a = np.zeros(0, dtype=np.complex128)               # Newton coeffs

    beta = norm(v)
    v = v / beta

    for s in range(maxrestart):

        arnoldi_vecs, Hess, Ritz, m = _arnoldi(
            A, dt, v, m_max, inner=inner, norm=norm)
        if m < m_max:
            logger.warn("Arnoldi only returned order %d instead of the "
                        "requested %d", m, m_max)
        if m == 0 and s == 0:
            # The input state appears to be an eigenstate
            eig_val = beta * Hess[0, 0]
            phase = func(eig_val)  # dt is absorbed in eig_val
            w = phase * v
            break

        # normalize Ritz points
        if s == 0:
            radius = _normalize_points(Ritz)
            center = 0.0
        assert(radius > 0.0), "Radius is zero"

        # get Leja points (i.e. Ritz points in the proper order)
        n_s = len(Z)
        Z, m = _extend_leja(Z, Ritz, m)  # Z now contains m new Leja points
        assert(m > 0), "No new Leja points"
        a = _extend_newton_coeffs(func, a, Z, center, radius)

        R = np.matrix(np.zeros(shape=(m+1, 1), dtype=np.complex128))
        R[0, 0] = beta
        P = a[n_s] * R
        for k in range(1, m):  # 1..m-1
            R = (np.dot(Hess, R) - Z[n_s+k-1] * R) / radius
            P += a[n_s+k] * R

        wp = zero(v)
        for i in range(m):  # 0 .. m-1
            wp += P[i, 0] * arnoldi_vecs[i]

        w += wp

        # starting vector for next iteration
        R = (np.dot(Hess, R) - Z[n_s+m-1] * R) / radius
        beta = np.linalg.norm(R)
        R /= beta
        # beta would be the norm of v, with the above normalization, v will now
        # be normalized
        v = zero(v)
        for i in range(m+1):  # 0 .. m
            v += R[i, 0] * arnoldi_vecs[i]

        if (beta*abs(a[-1])/(1+norm(w)) < tol):
            if logger.isEnabledFor(logging.DEBUG):  # pragma: nocover
                logger.debug("Converged at restart %s", s)
                logger.debug("norm of wp     : %s", norm(wp))
                logger.debug("norm of w      : %s", norm(w))
                logger.debug("beta           : %s", beta)
                logger.debug(
                    "|R*a[-1]|/|w|  : %s", np.linalg.norm(R) * a[-1] / norm(w))
                logger.debug("max Leja radius: %s", np.max(np.abs(Z)))
            break

        try:
            assert(not np.isnan(np.sum(v))), "v contains NaN"
            assert(not np.isnan(np.sum(w))), "w contains NaN"
        except (AttributeError, TypeError):
            pass

        if s == maxrestart - 1:
            logger.warn("DID NOT REACH CONVERGENCE")
            logger.warn("increase number of restarts")

    return w
