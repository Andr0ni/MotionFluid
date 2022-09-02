# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


import numpy as np
import time
from matplotlib import pyplot as plt

cur_time1 = 0.
x_max = 8
x_min = 0
n_x = 100
t_max = 1.2
m_t = 500
R00 = 1.3
g = 1.67
P0 = 50
h = (x_max - x_min) / n_x
t = t_max / m_t
p = np.zeros((n_x, m_t))
U = np.zeros((n_x, m_t))
R0 = np.zeros((n_x, m_t))

def start():
    for i in range(n_x):
        R0[i][0] = R00
        U[i][0] = 0
        p[i][0] = P(i * h, 4, 0.3)

    for i in range(m_t):
        U[0][i] = 0
        U[n_x - 1][i] = 0


def P(x, x0, r0):
    return P0 * np.exp(-((x - x0) ** 2) / r0 ** 2)


def draw(p0, u, ro):
    plt.ion()

    z = np.arange(x_min, x_max, h)
    y = p0.transpose()[0]
    fig, (p1, u1, ro1) = plt.subplots(3)
    fig.set_size_inches(10, 11)
    lineU, = u1.plot(z, y)
    lineRo, = ro1.plot(z, y)
    lineP, = p1.plot(z, y)
    fig.suptitle(f"gamma = {g}, Start density = {R00}, Start pressure = {P0}")


    for P, U, R in zip(p0.transpose(), u.transpose(), ro.transpose()):
        new_y = P

        p1.set_ylim(ymin=P.min(), ymax=P.max())
        lineP.set_xdata(np.arange(x_min, x_max, h))
        lineP.set_ydata(new_y)
        p1.set_title("Pressure")

        new_y = U

        u1.set_ylim(ymin=U.min(), ymax=U.max())
        lineU.set_xdata(np.arange(x_min, x_max, h))
        lineU.set_ydata(new_y)
        u1.set_title("Velocity")

        new_y = R

        ro1.set_ylim(ymin=R.min(), ymax=R.max())
        lineRo.set_xdata(np.arange(x_min, x_max, h))
        lineRo.set_ydata(new_y)
        ro1.set_title("Density")

        fig.canvas.draw()
        fig.canvas.flush_events()

        time.sleep(0.1)


def konserv_laks():
    start()
    for i in range(m_t - 1):
        for j in range(1, n_x - 1):
            R0[j][i + 1] = (R0[j + 1][i] + R0[j - 1][i]) / 2 - t * (
                        R0[j + 1][i] * U[j + 1][i] - R0[j - 1][i] * U[j - 1][i]) / (2 * h)
            p[j][i + 1] = ((p[j + 1][i] + p[j - 1][i]) / 2 - U[j][i] * t * (p[j + 1][i] - p[j - 1][i]) / (2 * h) -
                        g * p[j][i] * t * (U[j + 1][i] - U[j - 1][i]) / (2 * h))
            U[j][i + 1] = (((R0[j + 1][i] * U[j + 1][i] - R0[j - 1][i] * U[j - 1][i]) / 2 - t * (R0[j + 1][i] *
                        (U[j + 1][i] ** 2) + p[j + 1][i] - R0[j - 1][i] * ( U[j - 1][i] ** 2) - p[j - 1][i]) /
                        (2 * h)) / R0[j][i + 1])
    draw(p, U, R0)


def nekonserv_laks():
    start()
    for i in range(m_t - 1):
        for j in range(1, n_x - 1):
            R0[j][i + 1] = ((R0[j + 1][i] + R0[j - 1][i]) / 2 - U[j][i] * t * (R0[j + 1][i] - R0[j - 1][i]) / (2 * h) -
                            R0[j][i] * t * (U[j + 1][i] - U[j - 1][i]) / (2 * h))
            p[j][i + 1] = ((p[j + 1][i] + p[j - 1][i]) / 2 - U[j][i] * t * (p[j + 1][i] - p[j - 1][i]) / (2 * h) -
                            g * p[j][i] * t * (U[j + 1][i] - U[j - 1][i]) / (2 * h))
            U[j][i + 1] = ((U[j + 1][i] + U[j - 1][i]) / 2 - U[j][i] * t * (U[j + 1][i] - U[j - 1][i]) / (2 * h) -
                            t * (p[j + 1][i] - p[j - 1][i]) / (2 * h * R0[j][i]))
    draw(p, U, R0)


def Euler():
    start()
    for i in range(m_t - 1):
        for j in range(1, n_x - 1):

            if U[j][i] > 0:
                a = 0
            else:
                a = 1

            R0[j][i + 1] = R0[j][i] - (1 - a) * U[j][i] * t * (R0[j][i] - R0[j - 1][i]) / h - a * U[j][
                i] * t * (R0[j + 1][i] - R0[j][i]) / h - R0[j][i] * t * (U[j + 1][i] - U[j - 1][i]) / (2 * h)
            p[j][i + 1] = p[j][i] - (1 - a) * U[j][i] * t * (p[j][i] - p[j - 1][i]) / h - a * U[j][
                i] * t * (p[j + 1][i] - p[j][i]) / h - g * p[j][i] * t * (U[j + 1][i] - U[j - 1][i]) / (2 * h)
            U[j][i + 1] = U[j][i] - (1 - a) * U[j][i] * t * (U[j][i] - U[j - 1][i]) / h - a * U[j][
                i] * t * (U[j + 1][i] - U[j][i]) / h - t * (p[j + 1][i] - p[j - 1][i]) / (2 * h * R0[j][i])
    draw(p, U, R0)


#konserv_laks()
#nekonserv_laks()
Euler()
