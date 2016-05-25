import sympy as sym
import math
import numpy as num

t, omega, tau, p = sym.symbols('t omega tau p')

def countE(p):
    return ((2 / p) ** 2) * (sym.sin(p / 2) ** 2)


def tailor():
    return (((2 / p) ** 2) * (sym.sin(p / 2) ** 2)).subs(p, math.pi) + \
	    sym.diff(((2 / p) ** 2) * (sym.sin(p / 2)) ** 2, p).subs(p, math.pi) * \
            (p + math.pi) + \
            sym.diff(((2 / p) ** 2) * (sym.sin(p / 2)) ** 2, p, 2).subs(p, math.pi) * \
            ((p + math.pi) ** 2) / math.factorial(2) + \
            sym.diff(((2 / p) ** 2) * (sym.sin(p / 2)) ** 2, p, 3).subs(p, math.pi) * \
            ((p + math.pi) ** 3) / math.factorial(3) + \
            sym.diff(((2 / p) ** 2) * (sym.sin(p / 2)) ** 2, p, 4).subs(p, math.pi) * \
            ((p + math.pi) ** 4) / math.factorial(4)


if __name__ == "__main__":
    omega = 1.
    tau = 0.01
    s = num.linspace(0, math.pi, int(math.ceil(math.pi / tau)))

    with open("dotsE.dat", "w") as f:
        for val in s:
            f.write("" + str(val) + " " + str(countE(val)) + "\n")

    with open("tailor.dat", "w") as f:
        f.write(str(tailor()) + "\n")

