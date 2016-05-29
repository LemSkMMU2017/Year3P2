# -*- coding: utf-8 -*-
import sympy as sym

V, t, U, omega, tau = sym.symbols('V t U omega tau')
f = None

def ode_source_term(u):
    """
    Возвращает функцию источника ОДУ, равную u'' + omega**2*u.
    u --- символьная функция от t."""
    return sym.diff(u(t), t, t) + omega ** 2 * u(t)

def residual_discrete_eq(u):
    """
    Возвращает невязку разностного уравнения на заданной u.
    """
    R =  DtDt(u, tau) + omega ** 2 * u(t) - f
    return sym.simplify(R)

def residual_discrete_eq_step1(u):
    """
    Возвращает невязку разностного уравнения на первом шаге
    на заданной u.
    """
    R = u(tau) - (2 - tau ** 2 * omega ** 2) * U / 2 - tau * V - tau ** 2 * f.subs(t, 0) / 2
    return sym.simplify(R)

def DtDt(u, tau):
    """
	Возвращает вторую разностную производную от u.
	u --- символьная функция от t.
    """
    return (u(t + tau) - 2 * u(t) + u(t - tau)) / (2 * tau)

if __name__ == '__main__':
    linear = lambda t: V * t + U
    f = ode_source_term(linear)

    print "Погрешность:"
    print residual_discrete_eq(linear)

    print "Погрешность на 1 временном слое:"
    print residual_discrete_eq_step1(linear)