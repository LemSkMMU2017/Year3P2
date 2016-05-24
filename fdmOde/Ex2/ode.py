# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt


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
    u[1] = u[0] - 0.5 * tau ** 2 * omega ** 2 * u[0]
    for n in range(1, Nt):
        u[n+1] = 2 * u[n] - u[n-1] - tau ** 2 * omega ** 2 * u[n]
    return u, t


def u_exact(t, U, omega):
    return U * np.cos(omega * t)


def visualize(u, t, U, omega):
    plt.plot(t, u, 'r--o')
    t_fine = np.linspace(0, t[-1], 1001)  # мелкая сетка для точного решения
    u_e = u_exact(t_fine, U, omega)
    plt.hold('on')
    plt.plot(t_fine, u_e, 'b-')
    plt.legend([u'приближенное', u'точное'], loc='upper left')
    plt.xlabel('$t$')
    plt.ylabel('$u$')
    tau = t[1] - t[0]
    plt.title('$\\tau = $ %g' % tau)
    umin = 1.2*u.min();  umax = -umin
    plt.axis([t[0], t[-1], umin, umax])
    plt.savefig('tmp1.png')
    plt.savefig('tmp1.pdf')


