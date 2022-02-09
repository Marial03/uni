import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

warnings.filterwarnings('ignore', 'The iteration is not making good progress')


def euler(u_0, v_0, t_n, c, d, h):
    m = int(t_n / h)  # количество узлов
    y_u = np.zeros((m + 1,))
    y_v = np.zeros((m + 1,))
    t = np.linspace(0, t_n, m + 1)
    y_u[0] = u_0
    y_v[0] = v_0
    for i in range(m):
        if y_v[i - 1] >= 30:
            y_v[i] = c
            y_u[i] = y_u[i] + d
        y_v[i + 1] = y_v[i] + h * dv_dt(y_u[i], y_v[i])
        y_u[i + 1] = y_u[i] + h * du_dt(y_u[i], y_v[i])
    return t, y_v


def implicit_euler(u_0, v_0, t_n, c, d, h):
    m = int(t_n / h)  # количество узлов
    y_u = np.zeros((m + 1,))
    y_v = np.zeros((m + 1,))
    t = np.linspace(0, t_n, m + 1)
    y_u[0] = u_0
    y_v[0] = v_0

    def Phi_v(v, u_i, v_i):
        return v - v_i - h * dv_dt(u_i, v_i)

    def Phi_u(u, u_i, v_i):
        return u - u_i - h * du_dt(u_i, v_i)

    for i in range(m):
        if y_v[i - 1] >= 30:
            y_v[i] = c
            y_u[i] = y_u[i] + d
        y_v[i + 1] = optimize.fsolve(Phi_v, y_v[i], args=(y_u[i], y_v[i]))
        y_u[i + 1] = optimize.fsolve(Phi_u, y_u[i], args=(y_u[i], y_v[i]))
    return t, y_v


def runge_kutta(u_0, v_0, t_n, c, d, h):
    m = int(t_n / h)  # количество узлов
    y_u = np.zeros((m + 1,))
    y_v = np.zeros((m + 1,))
    t = np.linspace(0, t_n, m + 1)
    y_u[0] = u_0
    y_v[0] = v_0
    result = []
    T = t_n
    for i in range(m):
        if y_v[i - 1] >= 30:
            y_v[i] = c
            y_u[i] = y_u[i] + d
        k1_v = h * dv_dt(y_u[i], y_v[i])
        k1_u = h * du_dt(y_u[i], y_v[i])
        k2_v = h * dv_dt(y_u[i] + h / 2, y_v[i] + k1_v / 2)
        k2_u = h * du_dt(y_u[i] + h / 2, y_v[i] + k1_u / 2)
        k3_v = h * dv_dt(y_u[i] + h / 2, y_v[i] + k2_v / 2)
        k3_u = h * du_dt(y_u[i] + h / 2, y_v[i] + k2_u / 2)
        k4_v = h * dv_dt(y_u[i] + h, y_v[i] + k3_v)
        k4_u = h * du_dt(y_u[i] + h, y_v[i] + k3_u)
        y_v[i + 1] = y_v[i] + (1 / 6) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v)
        y_u[i + 1] = y_u[i] + (1 / 6) * (k1_u + 2 * k2_u + 2 * k3_u + k4_u)
    return t, y_v


def dv_dt(u, v):
    return 0.04 * (v ** 2) + 5 * v + 140 - u + 5


def du_dt(u, v):
    return a * (b * v - u)


# TS mode euler
a = 0.02
b = 0.2
c = -65
d = 6
v_0 = c
u_0 = b * v_0
t_approx, y_v_approx = euler(u_0, v_0, 20, c, d, h=0.05)
plt.subplot(2, 2, 1)
plt.title("Явный метод Эйлера")
plt.plot(t_approx, y_v_approx, linewidth=2)
plt.grid()

# PS mode euler
b = 0.25
c = -65
d = 6
v_0 = c
u_0 = b * v_0
t_approx, y_v_approx = euler(u_0, v_0, 20, c, d, h=0.05)
plt.subplot(2, 2, 2)
plt.plot(t_approx, y_v_approx, linewidth=2)
plt.grid()

# C mode euler
b = 0.2
c = -50
d = 2
v_0 = c
u_0 = b * v_0
t_approx, y_v_approx = euler(u_0, v_0, 20, c, d, h=0.05)
plt.subplot(2, 2, 3)
plt.plot(t_approx, y_v_approx, linewidth=2)
plt.grid()

# FS mode euler
a = 0.1
c = -65
v_0 = c
u_0 = b * v_0
t_approx, y_v_approx = euler(u_0, v_0, 20, c, d, h=0.05)
plt.subplot(2, 2, 4)
plt.plot(t_approx, y_v_approx, linewidth=2)
plt.grid()
plt.show()

# RUNGE_KUTTA
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------


# TS mode runge-kutta
a = 0.02
b = 0.2
c = -65
d = 6
v_0 = c
u_0 = b * v_0
t_approx, y_v_approx = runge_kutta(u_0, v_0, 100, c, d, h=0.05)
plt.subplot(2, 2, 1)
plt.title("Метод Рунге-Кутта")
plt.plot(t_approx, y_v_approx, linewidth=2)
plt.grid()

# PS mode runge-kutta
b = 0.25
c = -65
d = 6
v_0 = c
u_0 = b * v_0
t_approx, y_v_approx = runge_kutta(u_0, v_0, 20, c, d, h=0.05)
plt.subplot(2, 2, 2)
plt.plot(t_approx, y_v_approx, linewidth=2)
plt.grid()

# C mode runge-kutta
b = 0.2
c = -50
d = 2
v_0 = c
u_0 = b * v_0
t_approx, y_v_approx = runge_kutta(u_0, v_0, 20, c, d, h=0.05)
plt.subplot(2, 2, 3)
plt.plot(t_approx, y_v_approx, linewidth=2)
plt.grid()

# FS mode runge-kutta
a = 0.1
c = -65
v_0 = c
u_0 = b * v_0
t_approx, y_v_approx = runge_kutta(u_0, v_0, 20, c, d, h=0.05)
plt.subplot(2, 2, 4)
plt.plot(t_approx, y_v_approx, linewidth=2)
plt.grid()
plt.show()

# IMPLICIT EULER
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

# TS mode euler
a = 0.02
b = 0.2
c = -65
d = 6
v_0 = c
u_0 = b * v_0
t_approx, y_v_approx = implicit_euler(u_0, v_0, 20, c, d, h=0.05)
plt.subplot(2, 2, 1)
plt.title("Невный метод Эйлера")
plt.plot(t_approx, y_v_approx, linewidth=2)
plt.grid()

# PS mode euler
b = 0.25
c = -65
d = 6
v_0 = c
u_0 = b * v_0
t_approx, y_v_approx = implicit_euler(u_0, v_0, 20, c, d, h=0.05)
plt.subplot(2, 2, 2)
plt.plot(t_approx, y_v_approx, linewidth=2)
plt.grid()

# C mode euler
b = 0.2
c = -50
d = 2
v_0 = c
u_0 = b * v_0
t_approx, y_v_approx = implicit_euler(u_0, v_0, 20, c, d, h=0.05)
plt.subplot(2, 2, 3)
plt.plot(t_approx, y_v_approx, linewidth=2)
plt.grid()

# FS mode euler
a = 0.1
c = -65
v_0 = c
u_0 = b * v_0
t_approx, y_v_approx = implicit_euler(u_0, v_0, 20, c, d, h=0.05)
plt.subplot(2, 2, 4)
plt.plot(t_approx, y_v_approx, linewidth=2)
plt.grid()
plt.show()
