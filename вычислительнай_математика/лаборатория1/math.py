import numpy as np
import math
from numpy import linalg as LA
import matplotlib.pyplot as plt

def punkt_one(x, y):
    qs_coeff = qubic_spline_coeff(x, y)
    X = np.linspace(0,1,1000)
    Y = np.zeros(1000)
    for i in range(0, 1000):
        Y[i] = qubic_spline(X[i], qs_coeff, x)
    plt.scatter(x, y)
    plt.plot(x, y)
    plt.plot(X, Y)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.show()


def qubic_spline_coeff(x_nodes, y_nodes):
    # Инициализация массива сплайнов
    n = len(x_nodes)
    A = np.zeros((n, n))
    h = np.zeros((n - 1))

    for i in range(0, n - 1):
        h[i] = x_nodes[i + 1] - x_nodes[i]

    for i in range(1, n - 1):
        A[i][i] = 2 * (h[i] + h[i - 1])
        A[i][i + 1] = h[i]
        A[i][i - 1] = h[i - 1]
    A[0][0] = 1
    A[n - 1][n - 1] = 1


    B = np.zeros((n, 1))
    for i in range(1, n - 1):
        B[i][0] = 3 / h[i] * (y_nodes[i + 1] - y_nodes[i]) - 3 / h[i - 1] * (y_nodes[i] - y_nodes[i - 1])
    A_inv = LA.inv(A)
    C = A_inv.dot(B)

    a = np.zeros((n - 1, 1))
    b = np.zeros((n - 1, 1))
    c = np.zeros((n - 1, 1))
    d = np.zeros((n - 1, 1))

    for i in range(0, n - 1):
        a[i] = y_nodes[i]
        c[i] = C[i]
        # h1 = x_nodes[i + 1] - x_nodes[i]
        b[i] = (y_nodes[i + 1] - y_nodes[i]) / h[i] - h[i] * (C[i + 1] + 2 * C[i]) / 3
        d[i] = (C[i + 1] - C[i]) / 3 / h[i]

    con = np.column_stack([a, b, c, d])
    print(con)
    return con


def qubic_spline(xx, qs_coeff, x_nodes):
    for i in range(0, len(x_nodes)-1):
        if xx < x_nodes[i + 1]:
            s = qs_coeff[i][0] + qs_coeff[i][1] * (xx - x_nodes[i]) + qs_coeff[i][2] * pow(xx - x_nodes[i], 2) + qs_coeff[i][3] * pow(xx - x_nodes[i], 3)
            return s
    s = qs_coeff[-1][0] + qs_coeff[-1][1] * (xx - x_nodes[-2]) + qs_coeff[-1][2] * pow(xx - x_nodes[-2], 2) + qs_coeff[-1][3] * pow(xx - x_nodes[-2], 3)
    return s


def d_qubic_spline(xx, qs_coeff, x_nodes):
    for i in range(0, len(x_nodes)-1):
        if xx > x_nodes[i]:
            d_s = qs_coeff[i][1] + 2 * qs_coeff[i][2] * xx + 2 * qs_coeff[i][2] * x_nodes[i] + 3 * qs_coeff[i][3] * pow(
                xx, 2) - 6 * qs_coeff[i][3] * xx * x_nodes[i] + 3 * qs_coeff[i][3] * pow(x_nodes[i], 2)
            return d_s
    d_s = qs_coeff[-1][1] + 2 * qs_coeff[-1][2] * xx + 2 * qs_coeff[-1][2] * x_nodes[-2] + 3 * qs_coeff[-1][3] * pow(
        xx, 2) - 6 * qs_coeff[-1][3] * xx * x_nodes[-2] + 3 * qs_coeff[-1][3] * pow(x_nodes[-2], 2)
    return d_s


def bazis_polinom(i, x, x_nodes):
    l = 1
    for j in range(0, len(x_nodes)):
        if i != j:
            l = l * ((x - x_nodes[j]) / (x_nodes[i] - x_nodes[j]))
    return l


def interpoi_polinom(x, x_nodes, y_nodes):
    L = 0
    for i in range(0, len(x_nodes)):
        L = L + y_nodes[i] * bazis_polinom(i, x, x_nodes)
    return L


def generate(x, n):
    vectora = np.zeros((n, len(x)))
    for i in range(0, n):
        for j in range(0, len(x)):
            vectora[i][j] = x[j] + np.random.normal(0, 0.01)
    return vectora

def interpolanty_Lagranja_y(x_nodes, y_nodes, x, n):
    vectora = generate(y_nodes, n)
    all_y = []
    y = np.zeros(n)
    for i in range(0, n):
        for j in range(0, n):
            y[j] = interpoi_polinom(x[j], x_nodes, vectora[i])
        all_y.append(y.copy())
        plt.plot(x, y)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.show()
    return all_y

def punkt_three(x_nodes, y_nodes):
    n = 1000
    x = np.linspace(0, 1, n)

    vectora = generate(x_nodes, n)
    all_y = []
    y = np.zeros(n)
    for i in range(0, n):
        for j in range(0, n):
            y[j] = interpoi_polinom(x[j], vectora[i], y_nodes)
        all_y.append(y.copy())
        plt.plot(x, y)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.show()

    h_u = np.zeros(n)
    h_l = np.zeros(n)
    middle_Y = np.zeros(n)
    for j in range(0, n):
        sum_y = 0
        cur_y = np.zeros(n)
        for i in range(0, n):
            cur_y[i] += all_y[i][j]
            sum_y += all_y[i][j]
        middle_Y[j] = sum_y / n
        cur_y.sort()
        h_u[j] = cur_y[math.floor(n * 0.05)]
        h_l[j] = cur_y[math.floor(n * 0.95)]
    plt.plot(x, h_u)
    plt.plot(x, h_l)
    plt.scatter(x_nodes, y_nodes)
    plt.plot(x, middle_Y)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.show()


def punct_four(x_nodes, y_nodes):
    n = 1000
    x = np.linspace(0, 1, n)

    vectora = generate(y_nodes, n)
    all_y = []
    y = np.zeros(n)
    for i in range(0, n):
        for j in range(0, n):
            y[j] = interpoi_polinom(x[j], x_nodes, vectora[i])
        all_y.append(y.copy())
        plt.plot(x, y)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.show()

    h_u = np.zeros(n)
    h_l = np.zeros(n)
    middle_Y = np.zeros(n)
    for j in range(0, n):
        sum_y = 0
        cur_y = np.zeros(n)
        for i in range(0, n):
            cur_y[i] += all_y[i][j]
            sum_y += all_y[i][j]
        middle_Y[j] = sum_y / n
        cur_y.sort()
        h_u[j] = cur_y[math.floor(n * 0.05)]
        h_l[j] = cur_y[math.floor(n * 0.95)]
    plt.plot(x, h_u)
    plt.plot(x, h_l)
    plt.scatter(x_nodes, y_nodes)
    plt.plot(x, middle_Y)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.show()


def five_three(x_nodes, y_nodes):
    n = 100
    x = np.linspace(0, 1, n)
    vectora = generate(x_nodes, n)
    all_y = []
    y = np.zeros(n)

    for i in range(0, n):
        qs_coeff = qubic_spline_coeff(vectora[i], y_nodes)

        for j in range(0, n):
            y[j] = qubic_spline(x[j], qs_coeff, vectora[i])
        all_y.append(y.copy())
        plt.plot(x, y)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.show()

    h_u = np.zeros(n)
    h_l = np.zeros(n)
    middle_Y = np.zeros(n)
    for j in range(0, n):
        sum_y = 0
        cur_y = np.zeros(n)
        for i in range(0, n):
            cur_y[i] += all_y[i][j]
            sum_y += all_y[i][j]
        middle_Y[j] = sum_y / n
        cur_y.sort()
        h_u[j] = cur_y[math.floor(n * 0.05)]
        h_l[j] = cur_y[math.floor(n * 0.95)]
    plt.plot(x, h_u)
    plt.plot(x, h_l)
    plt.plot(x, middle_Y)
    plt.scatter(x_nodes, y_nodes)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.show()


def five_four(x_nodes, y_nodes):
    n = 100
    x = np.linspace(0, 1, n)
    vectora = generate(y_nodes, n)
    all_y = []
    y = np.zeros(n)

    for i in range(0, n):
        qs_coeff = qubic_spline_coeff(x_nodes, vectora[i])

        for j in range(0, n):
            y[j] = qubic_spline(x[j], qs_coeff, x_nodes)
        all_y.append(y.copy())
        plt.plot(x, y)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.show()

    h_u = np.zeros(n)
    h_l = np.zeros(n)
    middle_Y = np.zeros(n)
    for j in range(0, n):
        sum_y = 0
        cur_y = np.zeros(n)
        for i in range(0, n):
            cur_y[i] += all_y[i][j]
            sum_y += all_y[i][j]
        middle_Y[j] = sum_y / n
        cur_y.sort()
        h_u[j] = cur_y[math.floor(n * 0.05)]
        h_l[j] = cur_y[math.floor(n * 0.95)]
    plt.plot(x, h_u)
    plt.plot(x, h_l)
    plt.plot(x, middle_Y)
    plt.scatter(x_nodes, y_nodes)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True)
    plt.show()


x = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
y = [3.37, 3.95, 3.73, 3.59, 3.15, 3.15, 3.05, 3.86, 3.6, 3.7, 3.02]

punkt_one(x,y)
punkt_three(x, y)
punct_four(x,y)
five_three(x,y)
five_four(x, y)
