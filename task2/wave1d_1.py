# -*- coding: utf-8 -*-

import numpy as np

def solver(I, V, f, c, l, tau, gamma, T, user_action=None):
    K = int(round(T/tau))
    t = np.linspace(0, K*tau, K+1)
    dx = tau*c/float(gamma)
    N = int(round(l/dx))
    x = np.linspace(0, l, N+1)
    C2 = gamma**2                    # вспомогательная переменная
    if f is None or f == 0 :
        f = lambda x, t: 0
    if V is None or V == 0:
        V = lambda x: 0

    y   = np.zeros(N+1)   # Массив с решением на новом временном слое n+1
    y_1 = np.zeros(N+1)   # Решение на предыдущем слое n
    y_2 = np.zeros(N+1)   # Решение на слое n-1

    import time;  t0 = time.clock()  # для измерения процессорного времени

    # Задаем начальное условие
    for i in range(0,N+1):
        y_1[i] = I(x[i])

    if user_action is not None:
        user_action(y_1, x, t, 0)

    # Используем специальную формулу для расчета на первом
    # временном шаге с учетом du/dt = 0
    n = 0
    for i in range(1, N):
        y[i] = y_1[i] + tau*V(x[i]) + \
               0.5*C2*(y_1[i-1] - 2*y_1[i] + y_1[i+1]) + \
               0.5*tau**2*f(x[i], t[n])
    y[0] = 0;  y[N] = 0

    if user_action is not None:
        user_action(y, x, t, 1)

    y_2[:] = y_1;  y_1[:] = y

    for n in range(1, K):
        for i in range(1, N):
            y[i] = - y_2[i] + 2*y_1[i] + C2*(y_1[i-1] - 2*y_1[i] + y_1[i+1]) + tau**2*f(x[i], t[n])

        y[0] = 0; y[N] = 0 # Задаем граничные условия
        if user_action is not None:
            if user_action(y, x, t, n+1):
                break
        y_2[:] = y_1;  y_1[:] = y

    cpu_time = t0 - time.clock()
    return y, x, t, cpu_time

def test_quadratic():

    def u_exact(x, t):
        return x*(l-x)*(1 + 0.5*t)

    def I(x):
        return u_exact(x, 0)

    def V(x):
        return 0.5*u_exact(x, 0)

    def f(x, t):
        return 2*(1 + 0.5*t)*c**2

    l = 2.5
    c = 1.5
    gamma = 0.75
    N = 6
    tau = gamma*(l/N)/c
    T = 18

    def assert_no_error(u, x, t, n):
        u_e = u_exact(x, t[n])
        diff = np.abs(u - u_e).max()
        print "t = ", t[n], "   diff = ", diff;
        tol = 1E-13
        assert diff < tol

    solver(I, V, f, c, l, tau, gamma, T,
           user_action=assert_no_error)

def viz(
    I, V, f, c, l, tau, gamma, T,
    umin, umax,
    animate=True,
    tool='matplotlib',
    solver_function=solver,
    ):

    all_u = []

    def plot_u_st(u, x, t, n):
        plt.plot(x, u, 'r-',
                 xlabel='x', ylabel='u',
                 axis=[0, l, umin, umax],
                 title='t=%f' % t[n], show=True)
        all_u.append(u)
        time.sleep(0.01) if t[n] == 0 else time.sleep(0.01)
        plt.savefig('tmp_%04d.png' % n)

    class PlotMatplotlib:
        def __call__(self, u, x, t, n):
            if n == 0:
                plt.ion()
                self.lines = plt.plot(x, u, 'r-')
                plt.xlabel('x');  plt.ylabel('u')
                plt.axis([0, l, umin, umax])
                plt.legend(['t=%f' % t[n]], loc='lower left')
            else:
                self.lines[0].set_ydata(u)
                plt.legend(['t=%f' % t[n]], loc='lower left')
                plt.draw()
            time.sleep(2) if t[n] == 0 else time.sleep(0.2)
            plt.savefig('tmp_%04d.png' % n)

    if tool == 'matplotlib':
        import matplotlib.pyplot as plt
        plot_u = PlotMatplotlib()
    elif tool == 'scitools':
        import scitools.std as plt
        plot_u = plot_u_st
    import time, glob, os


    for filename in glob.glob('tmp_*.png'):
        os.remove(filename)

    user_action = plot_u if animate else None
    u, x, t, cpu = solver_function(
        I, V, f, c, l, tau, gamma, T, user_action)

    fps = 4
    codec2ext = dict(flv='flv', libx264='mp4', libvpx='webm',
                     libtheora='ogg')
    filespec = 'tmp_%04d.png'
    movie_program = 'ffmpeg'
    for codec in codec2ext:
        ext = codec2ext[codec]
        cmd = '%(movie_program)s -r %(fps)d -i %(filespec)s '\
              '-vcodec %(codec)s movie.%(ext)s' % vars()
        os.system(cmd)

    if tool == 'scitools':
        plt.movie('tmp_*.png', encoder='html', fps=fps,
                  output_file='movie.html')
    return cpu, all_u

def guitar(gamma):
    l = 0.75
    x0 = 0.8*l
    a = 0.005
    freq = 440
    wavelength = 2*l
    c = freq*wavelength
    omega = 2*np.pi*freq
    num_periods = 1
    T = 2*np.pi/omega*num_periods
    tau = l/50./c

    def I(x):
        return a*x/x0 if x < x0 else a/(l-x0)*(l-x)

    umin = -1.2*a;  umax = -umin
    cpu, all_u = viz(I, 0, 0, c, l, tau, gamma, T, umin, umax,
              animate=True, tool='scitools')
    print "all_u: "
    all_u = np.array(all_u).reshape(len(all_u), -1);
    print all_u

if __name__ == '__main__':
    test_quadratic()
    import sys
    try:
        gamma = float(sys.argv[1])
        print u'Число Куранта gamma=%g' % gamma
    except IndexError:
        gamma = 0.85
    print u'Число Куранта: %.2f' % gamma
    guitar(gamma)

