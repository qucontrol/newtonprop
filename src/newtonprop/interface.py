"""Recommended interface for a stateful Newton propagator

Typically, an efficient propagator will rely on low-level data structures.
Since conversion between high- and low-level structures in every call of
:func:`newtonprop.propagator.step` would be wasteful, it is better to use a
stateful propagator that only does the conversion once, holds the required
state, and converts back only on demand.

While a stateful propagator could be achieved through a wrapper function and
closures, a more straightfoward implementation holds the necessary state in a
class. To encourage a consistent interface, you should subclass
:class:`NewtonPropagatorBase` for this purpose.
"""
from abc import ABCMeta, abstractmethod

import numpy as np

from .propagator import _step, _exp, _expmi, _expi


__all__ = ["NewtonPropagatorBase"]


class NewtonPropagatorBase(metaclass=ABCMeta):
    """Stateful interface for using the Newton propagator

    This must be subclassed before it can be used.

    Args:
        data: data that underlies the abstract operator $A$. The propagator
            will use this together with a time parameter `t` to construct the
            argument `A` for :func:`newtonprop.propagator.step`.
        func (callable or str): The scalar version of the operator-function
            $f$.  Must take one complex argument and return a
            complex value. If given as a string, the corresponding option will
            be looked up in :attr:`funcs`.
        m_max (int): Maximal Krylov dimension
        maxrestart (int): Maximum number of Newton restarts
        tol (float): Desired precision of the result

    Subclasses will customize the behavior mostly by overriding private
    methods. See the source code of :class:`NewtonPropagatorBase` for details.
    """

    #: Registry of `func` names
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
        """The current state

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
        # This method may convert the initial state v into a more efficient
        # internal format
        return v

    def _convert_state_from_internal(self, v):  # pragma: nocover
        # If `self._v` uses an internal format different from the `v` in
        # `set_initial_state`, this method most be overridden to convert
        # `self._v` back to the original format.
        return v

    @staticmethod
    def _convert_data_to_internal(data):  # pragma: nocover
        return data

    @abstractmethod
    def _inner(self, a, b):  # pragma: nocover
        # Inner product between two states in the format of `self._v`
        pass

    @abstractmethod
    def _norm(self, a):  # pragma: nocover
        # Norm of a state in the format of `self._v`
        pass

    @abstractmethod
    def _zero(self, a):  # pragma: nocover
        # A new-zero state compatible with the format of `self._v`
        pass

    @abstractmethod
    def _A(self, t):  # pragma: nocover
        # return a callable `A(v)` that takes a state `v` in the
        # format of `self._v` and returns a state in the same format
        pass
