from numpy import matrix, array, diag, float64
from math import sqrt

def svd2(M):
    E = (M[0, 0] + M[1, 1]) * 0.5
    F = (M[0, 0] - M[1, 1]) * 0.5
    G = (M[0, 1] + M[1, 0]) * 0.5
    H = (M[0, 1] - M[1, 0]) * 0.5

    # one half of w1 + w2
    hw1pw2 = sqrt(E * E + H * H)
    # one half of w1 - w2
    hw1mw2 = sqrt(F * F + G * G)

    # cos(y - b)
    cymb = F / hw1mw2
    # cos(y + b)
    cypb = E / hw1pw2

    cc = sqrt((1 + cymb) * (1 + cypb))
    ss = sqrt((1 - cymb) * (1 - cypb))
    cs = sqrt((1 + cymb) * (1 - cypb))
    sc = sqrt((1 - cymb) * (1 + cypb))

    cb = (cc - ss) * 0.5
    sb = (sc + cs) * 0.5

    u = matrix([[cb, -sb], [sb, cb]], dtype=float64)
    d = array([hw1pw2 + hw1mw2, hw1pw2 - hw1mw2], dtype=float64)
    vT = matrix(diag(1.0 / d)).dot(u.T).dot(M)

    return u, d, vT

def pinv2(M):
    EPS = 1e-12
    u, s, vT = svd2(M)

    print(s)
    for k in range(s.size):
        if abs(s[k]) < EPS:
            s[k] = 0.0
        else:
            s[k] = 1.0 / s[k]

    return vT.T.dot(matrix(diag(s))).dot(u.T)