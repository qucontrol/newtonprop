"""Tests for `newtonprop` package."""
from pkg_resources import parse_version

from math import pi
import numpy as np
import pytest

import newtonprop

N = 10


def test_valid_version():
    """Check that the package defines a valid __version__"""
    assert parse_version(newtonprop.__version__) >= parse_version("0.1.0")


@pytest.fixture
def liouvillian():

    two_pi = 2.0 * pi
    w_c = two_pi * 0.01        # GHz
    gamma = two_pi * 0.1       # GHZ
    E0 = 0.1                   # GHZ

    # Drift Hamiltonian
    H0 = np.array(np.zeros(shape=(N, N), dtype=np.complex128))
    for i in range(N):
        H0[i, i] = i * w_c

    # Control Hamiltonian
    H1 = np.array(np.zeros(shape=(N, N), dtype=np.complex128))
    for i in range(1, N):
        H1[i-1, i] = np.sqrt(float(i))
        H1[i, i-1] = H1[i-1, i]

    # Total Hamiltonian
    H = H0 + E0 * H1

    # Dissipator
    a = np.array(np.zeros(shape=(N, N), dtype=np.complex128))
    for i in range(1, N):
        a[i-1, i] = np.sqrt(i)

    a_dag = a.conj().T

    def L(rho):
        Lrho = -1j * (H@rho - rho@H)
        Lrho += gamma * (
            a @ rho @ a_dag - 0.5 * (a_dag @ a @ rho + rho @ a_dag @ a))
        return Lrho

    return L


@pytest.fixture
def rho0():
    rho = np.matrix(np.zeros(shape=(N, N), dtype=np.complex128))
    rho[:, :] = 1.0 / float(N)
    return rho


def test_newton(liouvillian, rho0):
    assert abs(1-np.trace(rho0)) < 1e-12
    rho = newtonprop.propagator.step(liouvillian, rho0, dt=0.1, func=np.exp)
    assert abs(1-np.trace(rho)) < 1e-12
