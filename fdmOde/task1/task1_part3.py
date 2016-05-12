import sympy as sym

# define symbols variable
V, t, U, omega, tau = sym.symbols('V t U omega tau')

def ode_source_term(u):
    """calculate f function of
       Lu(t) = f(t)
    """

    return sym.diff(u(t), t, t) + omega ** 2 * u(t)

def residual_discrete_eq(u, f):
    """calculate deficiency of u(t) solution"""

    R =  DtDt(u, tau) + omega ** 2 * u(t) - f
    return sym.simplify(R)

def residual_discrete_eq_step1(u, f):
    """calculate deficiency of u(t) solution in the first step"""

    R = u(tau) - (2 - tau ** 2 * omega ** 2) * U / 2 - tau * V  - tau ** 2 * f.subs(t, 0) / 2
    return sym.simplify(R)

def DtDt(u, tau):
    """calculate u(t) on the difference scheme"""

    return (u(t + tau) - 2 * u(t) + u(t - tau)) / (2 * tau)

if __name__ == '__main__':
    linear = lambda t: V * t + U
    f = ode_source_term(linear)

    print "Check ODE deficiency:"
    print residual_discrete_eq(linear, f)

    print "Check ODE deficiency in the first step:"
    print residual_discrete_eq_step1(linear, f)
