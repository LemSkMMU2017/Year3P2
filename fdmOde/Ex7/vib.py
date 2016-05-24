# -*- coding: utf-8 -*-

import numpy as np
# import matplotlib.pyplot as plt
import scitools.std as plt


def solver(U, V, m, b, s, F, tau, T, damping='linear'):
    """
    Решает задачу m*u'' + f(u') + s(u) = F(t) for t in (0,T],
    u(0)=U и u'(0)=V,
    конечно-разностной схемой с шагом tau.
    Если затухание 'линейно', то f(u')=b*u, если затухание 'квадратичное', 
    то f(u')=b*u'*abs(u').
    F(t) и s(u) --- функции Python.
    """
    tau = float(tau);
    b = float(b);
    m = float(m)  # avoid integer div.
    N = int(round(T / tau))
    u = np.zeros(N + 1)
    t = np.linspace(0, N * tau, N + 1)

    uBot = U
    uTop = ''
    uCur = ''
    if damping == 'linear':
        uCur = uBot + tau * V + tau ** 2 / (2 * m) * (-b * V - s(uBot) + F(t[0]))
    elif damping == 'quadratic':
        uCur = uBot + tau * V + \
               tau ** 2 / (2 * m) * (-b * V * abs(V) - s(uBot) + F(t[0]))

    with open('funcU.txt', 'w') as f:
    	for n in range(1, N):
        	if damping == 'linear':
            		uTop = (2 * m * uCur + (b * tau / 2 - m) * uBot + \
                    tau ** 2 * (F(t[n]) - s(uCur))) / (m + b * tau / 2)
        	elif damping == 'quadratic':
            		uTop = (2 * m * uCur - m * uBot + b * uCur * abs(uCur - uBot)\
                    	+ tau ** 2 * (F(t[n]) - s(uCur))) / \
                   	(m + b * abs(uCur - uBot))
            	f.write(str(t[n]) + ' ' + str(uTop) + '\n')
        	uBot = uCur
        	uCur = uTop


def visualize(u, t, title='', filename='tmp'):
    plt.plot(t, u, 'b-')
    plt.xlabel('t')
    plt.ylabel('u')
    tau = t[1] - t[0]
    plt.title('tau=%g' % tau)
    umin = 1.2 * u.min();
    umax = 1.2 * u.max()
    plt.axis([t[0], t[-1], umin, umax])
    plt.title(title)
    plt.savefig(filename + '.png')
    plt.savefig(filename + '.pdf')
    plt.show()


import sympy as sym


def lhs_eq(t, m, b, s, u, damping='linear'):
    """Возвращает левую часть дифференциального уравнения как выражение sympy."""
    v = sym.diff(u, t)
    if damping == 'linear':
        return m * sym.diff(u, t, t) + b * v + s(u)
    else:
        return m * sym.diff(u, t, t) + b * v * sym.Abs(v) + s(u)


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--U', type=float, default=1.0)
    parser.add_argument('--V', type=float, default=0.0)
    parser.add_argument('--m', type=float, default=1.0)
    parser.add_argument('--b', type=float, default=0.0)
    parser.add_argument('--s', type=str, default='u')
    parser.add_argument('--F', type=str, default='0')
    parser.add_argument('--tau', type=float, default=0.05)
    parser.add_argument('--T', type=float, default=10)
    parser.add_argument('--window_width', type=float, default=30.,
                        help='Number of periods in a window')
    parser.add_argument('--damping', type=str, default='linear')
    parser.add_argument('--savefig', action='store_true')
    parser.add_argument('--SCITOOLS_easyviz_backend', default='matplotlib')
    a = parser.parse_args()
    from scitools.std import StringFunction
    s = StringFunction(a.s, independent_variable='u')
    F = StringFunction(a.F, independent_variable='t')
    U, V, m, b, tau, T, window_width, savefig, damping = \
        a.U, a.V, a.m, a.b, a.tau, a.T, a.window_width, a.savefig, \
        a.damping

    N = int(round(T / tau))
    u = np.zeros(N)
    t = np.zeros(N)
    solver(U, V, m, b, s, F, tau, T, damping)
    f = open('funcU.txt', 'r')
    i = 0
    for line in f:
        t[i], u[i] = line.split(' ')
        i += 1


def plot_empirical_freq_and_amplitude(u, t):
    minima, maxima = minmax(t, u)
    p = periods(maxima)
    a = amplitudes(minima, maxima)
    plt.figure()
    from math import pi
    omega = 2 * pi / p
    plt.plot(range(len(p)), omega, 'r-')
    plt.hold('on')
    plt.plot(range(len(a)), a, 'b-')
    ymax = 1.1 * max(omega.max(), a.max())
    ymin = 0.9 * min(omega.min(), a.min())
    plt.axis([0, max(len(p), len(a)), ymin, ymax])
    plt.legend(['estimated frequency', 'estimated amplitude'],
               loc='upper right')
    return len(maxima)


def visualize_front(u, t, window_width, savefig=False):
    """
	Визуализация приближенного и точного решений с использованием
	moving plot window и непрерывное изменение кривых от времени.
    P - приближенное значение периода.
    """
    import scitools.std as st
    from scitools.MovingPlotWindow import MovingPlotWindow

    umin = 1.2 * u.min();
    umax = -umin
    plot_manager = MovingPlotWindow(
            window_width=window_width,
            tau=t[1] - t[0],
            yaxis=[umin, umax],
            mode='continuous drawing')
    for n in range(1, len(u)):
        if plot_manager.plot(n):
            s = plot_manager.first_index_in_plot
            st.plot(t[s:n + 1], u[s:n + 1], 'r-1',
                    title='t=%6.3f' % t[n],
                    axis=plot_manager.axis(),
                    show=not savefig)  # drop window if savefig
            if savefig:
                print 't=%g' % t[n]
                st.savefig('tmp_vib%04d.png' % n)
    plot_manager.update(n)


def visualize_front_ascii(u, t, fps=10):
    """
    Визуализация приближенного и точного решений в коне 
    терминала (используются только символы ascii).
    """
    from scitools.avplotter import Plotter
    import time
    umin = 1.2 * u.min();
    umax = -umin

    p = Plotter(ymin=umin, ymax=umax, width=60, symbols='+o')
    for n in range(len(u)):
        print p.plot(t[n], u[n]), '%.2f' % (t[n])
        time.sleep(1 / float(fps))


def minmax(t, u):
    """
    Вычисляем все локальные минимумы и максимумы приближенного
    решения u(t), заданного массивами  u и t.
    Возвращает список минимумов и максимумов вида (t[i],u[i]).
    """
    minima = [];
    maxima = []
    for n in range(1, len(u) - 1, 1):
        if u[n - 1] > u[n] < u[n + 1]:
            minima.append((t[n], u[n]))
        if u[n - 1] < u[n] > u[n + 1]:
            maxima.append((t[n], u[n]))
    return minima, maxima


def periods(extrema):
    """
    По заданному списку точек минимума и максимума (t,u),
    возвращает массив соотвествующих локльных периодов.
    """
    p = [extrema[n][0] - extrema[n - 1][0]
         for n in range(1, len(extrema))]
    return np.array(p)


def amplitudes(minima, maxima):
    """
    По заданным спискам точек минимума и максимума (t,u), 
    возвращает массив соответствующих локальных амплитуд.
    """
    # Сравниваем первый максимум с первым минимумом и т.д.
    a = [(abs(maxima[n][1] - minima[n][1])) / 2.0
         for n in range(min(len(minima), len(maxima)))]
    return np.array(a)


if __name__ == '__main__':
    main()
    # test_constant()
    # test_sinusoidal()
    # test_mms()
    # test_quadratic()
