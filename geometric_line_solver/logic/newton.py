import numpy as np

EPS = 1e-9
MAX_ITER = 10000
ERROR_MAX_ITER = "Превышено число итераций метода Ньютона"


def newtons_method(get_jf, start_v):
    # число итераций
    k = 0
    cur_v = start_v
    while k < MAX_ITER:
        # расчет F и J
        j_matrix, f_vector = get_jf(cur_v)
        # критерий останова по F(V)
        # S = np.max(np.abs(delta_vector))
        # print('S = {}'.format(S))
        # if S <= EPS:
        #     return cur_v
        # j_matrix += np.diag(np.ones(j_matrix.shape[0], dtype=np.double) * 1e-18, )

        # решение СЛАУ
        try:
            delta_vector = np.linalg.solve(j_matrix, -f_vector)
        except np.linalg.LinAlgError as e:
            raise

        # прибавление поправок к текущим координатам V_k+1 = V_k + dV_k
        cur_v += delta_vector

        S = np.sqrt(np.sum(delta_vector ** 2))
        print('S = {}'.format(S))
        if S <= EPS:
            print('--------------------------------')
            return cur_v
        k = k + 1
    raise RuntimeError(ERROR_MAX_ITER)

# тесты
# def get_sin_for_newton(x):
#     return np.array([-np.sin(x)]).reshape(-1, 1), np.cos(x)
#
# def get_test_f_for_newton(x):
#     return np.array([1-np.cos(x)]).reshape(-1, 1), -(x - np.sin(x) - 0.25)