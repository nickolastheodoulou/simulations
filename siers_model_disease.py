#%%

from vpython import *

# GlowScript 2.9 VPython
# https://www.idmod.org/docs/hiv/model-seir.html
# https://en.wikipedia.org/wiki/Semi-implicit_Euler_method
# https://sci-hub.tw/10.1007/s10483-007-0914-x

# Define parameters.
b = 0.000001  # birth rate
mu = 0.0000005  # death rate
beta = 0.8  # rate of trasmission
alpha = 0.9  # the rate of disease infectivity
gamma = 0.90  # the rate of cure of disease
sigma = 0.01  # the rate of returning individuals vulnerable
nu = 0  # the ratio of the number of individuals who received the vaccine
q = 0  # giving treatment

def derivs(t, population):
    S = population[0]  # Fraction susceptible: Not yet infected and might become infected.
    I = population[1]  # Fraction infected: Ongoing disease.
    R = population[2]  # Fraction recovered: Gained immunity, isolated, or dead.
    E = population[3]
    N = S + I + R + E

    Sdot = (b * N) + (sigma * R) - (beta * S * I/ N) - (nu * S) - (mu * S)
    Edot = (beta * S * I/ N) - (alpha * E) - (mu * E)
    Idot = (alpha * E) + (q * I) - (gamma * I) - (mu * I)
    Rdot = (gamma * I) + (nu * S) - (sigma * R) - (mu * R)

    return [Sdot, Idot, Rdot, Edot]


# Set initial conditions. These are fractions, not whole numbers.
S = 50000000  # number susceptible: Not yet infected.
I = 500000  # number infected: Ongoing disease.
R = 1000  # number recovered: Gained immunity, isolated, or dead.
E = 5000000  # number exposed

# Create graphs.
s_plot = gcurve(color=color.blue, label="% susceptible")
i_plot = gcurve(color=color.red, label="% infected")
r_plot = gcurve(color=color.green, label="% recovered")
e_plot = gcurve(color=color.orange, label="% exposed")

time = 0
tmax = 1000
dt = 0.1


def mult(l1, scal):
    out = []
    for i in range(0, len(l1)):
        out.append(l1[i] * scal)
    return out


def add(l1, l2):
    out = []
    for i in range(0, len(l1)):
        out.append(l1[i] + l2[i])
    return out


population = [S, I, R, E]


while time < tmax:
    rate(50)
    # Use RK4 method.
    k1 = derivs(time, population)
    k2 = derivs(time + dt / 2, add(population, mult(k1, 0.5)))
    k3 = derivs(time + dt / 2, add(population, mult(k2, 0.5)))
    k4 = derivs(time + dt, add(population, k3))
    population = add(population,
                     mult(add(mult(k1, 1 / 6), add(mult(k2, 1 / 3), add(mult(k3, 1 / 3), mult(k4, 1 / 6)))), dt))
    s_plot.plot(time, population[0])
    i_plot.plot(time, population[1])
    r_plot.plot(time, population[2])
    e_plot.plot(time, population[3])

    time += dt
