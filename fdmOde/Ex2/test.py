import numpy as np
import ode
import unittest


TAU_SEQUENCE = [0.1, 0.5, 0.01]
OMEGA_SEQUENCE = [0.7, 0.5]
U_SEQUENCE = [3., 2., 1.]
T = 1


class Test(unittest.TestCase):
    def test_solution(self):
        tol = 1E-13
        for tau in TAU_SEQUENCE:
            for omega in OMEGA_SEQUENCE:
                for U in U_SEQUENCE:
                    u, t = ode.solver(U, omega, tau, T)
                    u_e = Test.__exact_numerical_u__(U, omega, t)
                    error_norm = Test.__error_norm__(u_e, u, tau)
                    self.assertLessEqual(error_norm, tol)

    @staticmethod
    def __exact_numerical_u__(U, omega, t):
        tau = t[1] - t[0]
        approximated_omega = Test.__approximate_omega__(omega, tau)
        return np.array([U * np.cos(approximated_omega * t_n) for t_n in t])

    @staticmethod
    def __approximate_omega__(omega, tau):
        return 2. * np.arcsin(omega * tau / 2.) / tau

    @staticmethod
    def __error_norm__(u_e, u, tau):
        return np.sqrt(tau * np.sum((u_e-u) ** 2))

if __name__ == '__main__':
    unittest.main()
