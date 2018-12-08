"""Recommended interface for a stateful Newton propagator

Typically, an efficient propagator will rely on low-level data structures.
Since conversion between high- and low-level structures in every call of
:func:`newtonprop.propagator.step` would be wasteful, it is better to use a
stateful propagator that only does the conversion once, holds the required
state, and converts back only on demand.

While a stateful propagator could be achieved through a wrapper function and
closures, a more straightforward implementation holds the necessary state in a
class. To encourage a consistent interface, you should subclass
:class:`NewtonPropagatorBase` for this purpose.
"""
from abc import ABCMeta, abstractmethod

import numpy as np

from .propagator import _step, _exp, _expmi, _expi


__all__ = ["NewtonPropagatorBase"]


class NewtonPropagatorBase(metaclass=ABCMeta):
    """Abstract base class for a stateful interface for using the Newton
    propagator

    This must be sub-classed before it can be used.

    Args:
        data: data that underlies the abstract operator $A$. The propagator
            will use this together with a time parameter `t` to construct the
            argument `A` for :func:`newtonprop.propagator.step`.
        func (callable or str): The scalar version of the operator-function
            $f$.  Must take one complex argument and return a
            complex value. If given as a string, the corresponding callable
            will be looked up in :attr:`funcs`.
        m_max (int): Maximal Krylov dimension. Initializes the corresponding
            attribute.
        maxrestart (int): Maximum number of Newton restarts. Initializes the
            corresponding attribute.
        tol (float): Desired precision of the result. Initializes the
            corresponding attribute.

    Subclasses will customize the behavior mostly by overriding a number of
    private methods, outlined below. These methods may act on the following
    (private) attributes:

    .. py:attribute:: _data

        Internal data that may be used to calculate $A$. Obtained from
        `data` via :meth:`_convert_data_to_internal`.

    .. py:attribute:: _v0

        The state set by :meth:`set_initial_state`, in its original format

    .. py:attribute:: _v

        The current state, in an internal format. The internal format is
        determined by :meth:`_convert_state_to_internal`.

    **Private methods:**

    .. automethod:: _check

    .. automethod:: _convert_state_to_internal

    .. automethod:: _convert_state_from_internal

    .. automethod:: _convert_data_to_internal

    .. automethod:: _inner

    .. automethod:: _norm

    .. automethod:: _zero

    .. automethod:: _A

    **Public methods, attributes, and properties:**

    .. py:attribute:: func

        The scalar version of the operator-function $f$, as a callable

    .. py:attribute:: m_max

        Maximal Krylov dimension (int)

    .. py:attribute:: maxrestart

        Maximum number of Newton restarts (int)

    .. py:attribute:: tol

        Desired precision of the result (float)
    """

    #: Registry of `func` names (class attribute)
    #:
    #: Cf. the `func` argument of :func:`newtonprop.propagator.step`. You may
    #: extend this with your own function names.
    funcs = {"exp": _exp, "expmi": _expmi, "expi": _expi}

    def __init__(self, data, func="exp", m_max=10, maxrestart=100, tol=1e-12):
        if isinstance(func, str):
            func = self.funcs[func]
        self._data = self._convert_data_to_internal(data)
        self.func = func
        self._v0 = None
        self._v = None
        self.m_max = m_max
        self.maxrestart = maxrestart
        self.tol = tol

    def set_initial_state(self, v):
        """Set the initial state for the propagation."""
        # Note: `self._v` does not need to be the same as the input `v`. In
        # general, you want to convert `v` to a more efficient internal format
        # and store that.
        self._v0 = v
        self._v = self._convert_state_to_internal(v)
        self._check()

    @property
    def state(self):
        """The current state (read-only)

        After a call to :meth:`set_initial_state`, this is the state set by
        that call. After a call to :meth:`step`, it is the result of that
        propagation step.
        """
        res = self._convert_state_from_internal(self._v)
        assert isinstance(res, self._v0.__class__)
        return res

    def step(self, t, dt):
        r"""Perform a single propagation step

        Args:
            t (float): time value at which to evaluate the abstract operator
                $A(t)$
            dt (float): length of the time step

        Construct an abstract operator $A(t)$ based on the `data` the
        propagator was initialized with, and calculate

        .. math::

            v_\text{out} = f(A(t)\,dt) \, v_\text{in}

        where $v_\text{in}$ is the current :attr:`state`, and $f$ is the
        function specified by `func` during initialization.

        The state $v_\text{out}$ becomes $v_\text{in}$ for the next call to
        :meth:`step`. The time parameter `t` should generally be the midpoint
        of the interval covered by $dt$.

        The results of the propagation step may be read from :attr:`state`. The
        :meth:`step` method does not return anything.
        """
        A = self._A(t)
        self._v = _step(
            A,
            self._v,
            dt,
            self.func,
            m_max=self.m_max,
            maxrestart=self.maxrestart,
            tol=self.tol,
            inner=self._inner,
            norm=self._norm,
            zero=self._zero,
        )

    # subclasses should customize the behavior primarily by overriding the
    # following private methods:

    def _check(self):  # pragma: nocover
        """Check completeness and consistency of attributes

        After :meth:`set_initial_state`, check that all attributes are defined
        and are consistent. If not, raise an :exc:`AttributeError`.
        """
        try:
            N = np.prod(self._v.shape)
            if self.m_max >= N:
                raise ValueError(
                    "m_max must be smaller than the system dimension"
                )
        except AttributeError:
            pass

    def _convert_state_to_internal(self, v):  # pragma: nocover
        """Convert the initial state `v` into a more efficient internal format.

        This is used when setting :attr:`_v` from the `v` given in
        :meth:`set_initial_state`. By default, the state is passed on
        unchanged.
        """
        return v

    def _convert_state_from_internal(self, v):  # pragma: nocover
        """Convert `v` in internal format of :attr:`_v` to original format

        If :attr:`_v` uses an internal format different from the `v` in
        :meth:`set_initial_state`, this method must be overridden to convert
        :attr:`_v` back to the original format. It is used by :attr:`state` for
        this purpose. By default, the internal state is passed on unchanged.
        """
        return v

    @staticmethod
    def _convert_data_to_internal(data):  # pragma: nocover
        """Convert `data` to :attr:`_data`

        Optional conversion of the `data` passed to the initializer to a more
        efficient internal format.  By default, the `data` is stored unchanged.
        """
        return data

    @abstractmethod
    def _inner(self, a, b):  # pragma: nocover
        """Inner product between two states `a`, `b` in the format of
        :attr:`_v`

        Must return a complex number. This is an abstract method, and thus must
        be implemented by any subclass of :class:`NewtonPropagatorBase`.

        An example implementation is::

            def _inner(self, a, b):
                return numpy.vdot(a, b)
        """
        pass

    @abstractmethod
    def _norm(self, a):  # pragma: nocover
        """Norm of a state `a` in the format of :attr:`_v`

        Must return a complex number. This is an abstract method, and thus must
        be implemented by any subclass of :class:`NewtonPropagatorBase`.

        An example implementation is::

            def _norm(self, a):
                return numpy.linalg.norm(a)
        """
        pass

    @abstractmethod
    def _zero(self, a):  # pragma: nocover
        """A new zero-state compatible with the format of :attr:`_v`

        Must allocate and return a new state in the same internal format as
        :attr:`_v`, but containing zeros, so that the state can set in-place.
        This is an abstract method, and thus must be implemented by any
        subclass of :class:`NewtonPropagatorBase`.

        An example implementation is::

            def _zero(self, a):
                return numpy.zeros(shape=a.shape, dtype=a.dtype)
        """
        pass

    @abstractmethod
    def _A(self, t):  # pragma: nocover
        """Callable ``A(v)`` encoding the abstract operator $A(t)$.

        Must return a callable ``A(v)`` that takes a state `v` in the
        internal format of :attr:`_v` and returns a state in the same format.
        This is an abstract method, and thus must be implemented by any
        subclass of :class:`NewtonPropagatorBase`.

        Assuming ``self._data`` contains a nested list ``[L0, [L1, eps]]``
        where ``L0`` and ``L1`` are numpy matrices containing a superoperator,
        and ``eps`` is a callable with argument `t` for a time-dependent
        control, an example implementation is::

            def _A(self, t):
                Lmatrix = (
                    self._data[0] + self._data[1][0] + self._data[1][1](t))

                def L(v):
                    return Lmatrix @ v

                return L

        .. Note::

            This method does not calculate ``A(v)`` directly, but returns a
            callable that does. This is to take into account time-dependent
            operators, which :func:`newtonprop.propagator.step` does not know
            anything about.
        """
        pass
