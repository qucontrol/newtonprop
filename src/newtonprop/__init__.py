"""Pure Python reference implementation of the Newton propagator for quantum
dynamics."""

__version__ = '0.1.0'

from .propagator import step as newton

__all__ = ['newton']
