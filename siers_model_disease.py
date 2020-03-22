#%%

from vpython import *

# GlowScript 2.9 VPython
# https://www.idmod.org/docs/hiv/model-seir.html
# https://en.wikipedia.org/wiki/Semi-implicit_Euler_method
# https://sci-hub.tw/10.1007/s10483-007-0914-x

#%%

# Define parameters.
beta = 0.35  # Rate of disease transmission. Includes washing your hands!
gamma = 0.1  # Rate of loss of recovery. 0 --> Once in recovery, stay in recovery.
xi = 0.005  # Rate at which individuals return to susceptible state due to loss of immunity
sigma = 0.01  #rate of latent individuals becoming infectious

def derivs(t, population):
    S = population[0]  # Fraction susceptible: Not yet infected and might become infected.
    I = population[1]  # Fraction infected: Ongoing disease.
    R = population[2]  # Fraction recovered: Gained immunity, isolated, or dead.
    E = population[3]
    N = S + I + R + E

    Sdot = -(beta * S * I/N) + (xi * R)
    Edot = (beta * S * I/N) - (sigma * E)
    Idot = (sigma * E) - (gamma * I)
    Rdot = (gamma * I) - (xi * R)

    return [Sdot, Idot, Rdot, Edot]


# Set initial conditions. These are fractions, not whole numbers.
S = 2000  # number susceptible: Not yet infected.
I = 40  # number infected: Ongoing disease.
R = 0  # number recovered: Gained immunity, isolated, or dead.
E = 4

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
