# -*- coding: utf-8 -*-

import numpy as np
from math import pi


def solver(U, omega, tau, T):
    """
    Решается задача
    u'' + omega**2*u = 0 для t из (0,T], u(0)=U и u'(0)=0,
    конечно-разностным методом с постоянным шагом tau
    """
    tau = float(tau)
    Nt = int(round(T/tau))
    u = np.zeros(Nt+1)
    t = np.linspace(0, Nt*tau, Nt+1)

    u[0] = U
    u[1] = u[0] - 0.5*tau**2*omega**2*u[0]
    for n in range(1, Nt):
        u[n+1] = 2*u[n] - u[n-1] - tau**2*omega**2*u[n]
    return u, t


def u_exact(t, U, omega):
    return U*np.cos(omega*t)


def minmax(t, u):
    """
    Вычисляем все локальные минимумы и максимумы приближенного
    решения u(t), заданного массивами  u и t.
    Возвращает список минимумов и максимумов вида (t[i],u[i]).
    """
    minima = []; maxima = []
    for n in range(1, len(u)-1, 1):
        if u[n-1] > u[n] < u[n+1]:
            minima.append((t[n], u[n]))
        if u[n-1] < u[n] > u[n+1]:
            maxima.append((t[n], u[n]))
    return minima, maxima


def get_numerical_spades(extrema, num_periods):
    """
    По заданному списку точек максимума (t,u),
    возвращает массив точек, соответствующих пикам.
    """
    num = min(len(extrema), num_periods - 1)
    p = [extrema[n][0] for n in xrange(0, num)]
    return np.array(p)


def get_exact_spades(omega,  num_periods):
    P = 2 * pi / omega
    return np.array([P * n for n in xrange(1, num_periods)])


def get_deviations(exact_spades, numerical_spades):
    return exact_spades - numerical_spades


def get_ratio(deviations):
    return [deviations[m - 1] / float(m) for m in xrange(1, len(deviations))]


def test_linearity(U, omega, tau, num_periods):
    P = 2*pi/omega  # один период
    T = P*num_periods
    u, t = solver(U, omega, tau, T)
    minima, maxima = minmax(t, u)
    exact_spades = get_exact_spades(omega, num_periods)
    numerical_spades = get_numerical_spades(maxima, num_periods)
    deviations = get_deviations(exact_spades, numerical_spades)
    ratios = get_ratio(deviations)
    assert abs(ratios[-2] - ratios[-1]) <= tau / len(ratios)
    print 'OK/Erorr is linearly dependent on m/---U={}---omega={}---tau={}---num_periods={}---'.format(U, omega, tau, num_periods)

if __name__ == '__main__':
    test_linearity(1., 2 * pi, 0.05, 100)
    test_linearity(1., pi, 0.01, 100)
    test_linearity(3., 2 * pi, 0.001, 1000)

