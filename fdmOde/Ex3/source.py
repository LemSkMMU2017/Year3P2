import sympy as sym

# define symbols variable
# m u'' + f(u') + s(u) = F(t), u(0) = U, u'(0) = V, t in (0, T]
U, V, T, t, tau, m, b, c = sym.symbols('U V T t tau m b c')

def solver(U, V, m, b, s, F, tau, T, damping='linear'):
    tau = float(tau); b = float(b); m = float(m) # avoid integer div.
    N = int(round(T/tau))
    u = np.zeros(N+1)
    t = np.linspace(0, N*tau, N+1)

    u[0] = U
    if damping == 'linear':
        u[1] = u[0] + tau*V + tau**2/(2*m)*(-b*V - s(u[0]) + F(t[0]))
    elif damping == 'quadratic':
        u[1] = u[0] + tau*V + \
               tau**2/(2*m)*(-b*V*abs(V) - s(u[0]) + F(t[0]))

    for n in range(1, N):
        if damping == 'linear':
            u[n+1] = (2*m*u[n] + (b*tau/2 - m)*u[n-1] +
                      tau**2*(F(t[n]) - s(u[n])))/(m + b*tau/2)
        elif damping == 'quadratic':
            u[n+1] = (2*m*u[n] - m*u[n-1] + b*u[n]*abs(u[n] - u[n-1])
                      + tau**2*(F(t[n]) - s(u[n])))/\
                      (m + b*abs(u[n] - u[n-1]))
    return u, t

#def linear_oscillation(u, b):
#    """Return f(u') = b * u'  -  linear oscillation func."""
#
#    return sym.diff(u(t), t) * b
#
def ode_source_term(u, m, f, s):
    """Calculate and return F(t) from equation: mu'' + f(u') + s(u) = F(t)"""

    R = m * sym.diff(u(t), t, t) + f(sym.diff(u(t), t)) + s(u(t))
    return sym.simplify(R)

def residual_discrete_eq(u, m, f, s):
    """calculate deficiency of u(t) solution"""

    R =  m * DtDt(u, tau) + f(Dt(u, tau)) + s(u(t))
    return sym.simplify(R)

def residual_discrete_eq_for_square(u, m, b, s):
    """calculate deficiency of u(t) solution"""

    R =  m * DtDt(u, tau) + b * Dt_left(u, tau) * abs(Dt_right(u, tau)) + s(u(t))
    return sym.simplify(R)

def residual_discrete_eq_step1_for_linear(u, s, F):
    """calculate deficiency of u(t) solution in the first step"""

    #R = u(tau) - (2 - tau ** 2 * omega ** 2) * U / 2 - tau * V  - tau ** 2 * f.subs(t, 0) / 2
    R = u(tau) - u(0) - tau * V - tau ** 2 / (2 * m) * (F.subs(t, 0) - s(U) - b * V)
    return sym.simplify(R)

def residual_discrete_eq_step1_for_square(u, s, F):
    """calculate deficiency of u(t) solution in the first step"""

    #R = u(tau) - (2 - tau ** 2 * omega ** 2) * U / 2 - tau * V  - tau ** 2 * f.subs(t, 0) / 2
    R = u(tau) - u(0) - tau * V - tau ** 2 / (2 * m) * (F.subs(t, 0) - s(U) - b * V * abs(V))
    return sym.simplify(R)

def DtDt(u, tau):
    """calculate u(t) on the difference scheme"""

    return (u(t + tau) - 2 * u(t) + u(t - tau)) / (tau ** 2)

def Dt(u, tau):
    """calculate u'(t) on the difference scheme"""

    return (u(t + tau) -  u(t - tau)) / (2 * tau)

def Dt_left(u, tau):
    """calculate u'(t) on the difference scheme"""

    return (u(t + tau) - u(t)) / tau

def Dt_right(u, tau):
    """calculate u'(t) on the difference scheme"""

    return (u(t) - u(t - tau)) / tau

if __name__ == '__main__':
# m * u'' + f(u') + s(u) = F(t), u(0) = U, u'(0) = V, t in (0, T]
    f = lambda t: b * t
    s = lambda t: c * t

    linear = lambda t: V * t + U
    square = lambda t: b * t * t + V * t + U
    F = ode_source_term(linear, m, f, s)

    print "Check ODE deficiency with linear vibration for linear func:"
    print sym.simplify(residual_discrete_eq(linear, m, f, s) - F)

    print "Check ODE deficiency in the first step with linear vibration for linear func:"
    print residual_discrete_eq_step1_for_linear(linear, s, F)

    F = ode_source_term(square, m, f, s)
    print "Check ODE deficiency with linear vibration for square func:"
    print sym.simplify(residual_discrete_eq(square, m, f, s) - F)

    print "Check ODE deficiency in the first step with linear vibration for square func :"
    print residual_discrete_eq_step1_for_linear(square, s, F)

    f = lambda t: b * t * abs(t)
    F = ode_source_term(linear, m, f, s)
    print "Check ODE deficiency with square vibration for linear func:"
    print sym.simplify(residual_discrete_eq_for_square(linear, m, b, s) - F)

    print "Check ODE deficiency in the first step with square vibration for linear func:"
    print residual_discrete_eq_step1_for_square(linear, s, F)

    F = ode_source_term(square, m, f, s)
    print "Check ODE deficiency with square vibration for square func:"
    print sym.simplify(residual_discrete_eq_for_square(square, m, b, s) - F)

    print "Check ODE deficiency in the first step with square vibration for square func:"
    print residual_discrete_eq_step1_for_square(square, s, F)
