import vpython

# GlowScript 2.9 VPython
# https://www.idmod.org/docs/hiv/model-seir.html
# https://en.wikipedia.org/wiki/Semi-implicit_Euler_method
# https://sci-hub.tw/10.1007/s10483-007-0914-x


def derivs(t, population):

    # Set initial conditions.

    # parameters
    birth_rate = 0.000001  # birth rate
    death_rate = 0.0000005  # death rate
    rate_of_transmission = 0.8  # rate of transmission
    rate_of_incubation = 1.0  # the rate of latent individuals becoming infectious (average duration of incubation is 1/\sigma)
    rate_of_recovery = 0.90  # the rate of cure of disease
    rate_of_loss_of_immunity = 0.01  # the rate at which recovered individuals return to the susceptible statue due to loss of immunity.
    rate_of_vaccination = 0  # the ratio of the number of individuals who received the vaccine
    rate_of_treatment = 0  # giving treatment



    S = population[0]  # number susceptible: Not yet infected.
    I = population[1]  # number infected: Ongoing disease.
    R = population[2]  # number recovered: Gained immunity, isolated, or dead.
    E = population[3]  # number exposed
    total_population = S + I + R + E

    Sdot = (birth_rate * total_population) + (rate_of_loss_of_immunity * R) - (rate_of_transmission * S * I / total_population) - (rate_of_vaccination * S) - (death_rate * S)
    Edot = (rate_of_transmission * S * I / total_population) - (rate_of_incubation * E) - (death_rate * E)
    Idot = (rate_of_incubation * E) + (rate_of_treatment * I) - (rate_of_recovery * I) - (death_rate * I)
    Rdot = (rate_of_recovery * I) + (rate_of_vaccination * S) - (rate_of_loss_of_immunity * R) - (death_rate * R)

    return [Sdot, Idot, Rdot, Edot]


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


def main():

    pop_susceptible = 50000000  # number susceptible: Not yet infected.
    pop_infectious = 500000  # number infected: Ongoing disease.
    pop_recovered = 1000  # number recovered: Gained immunity, isolated, or dead.
    pop_exposed = 5000000  # number exposed
    vector_population = [pop_susceptible, pop_infectious, pop_recovered, pop_exposed]

    # Create graphs.
    s_plot = vpython.gcurve(color=vpython.color.blue, label="% susceptible")
    i_plot = vpython.gcurve(color=vpython.color.red, label="% infected")
    r_plot = vpython.gcurve(color=vpython.color.green, label="% recovered")
    e_plot = vpython.gcurve(color=vpython.color.orange, label="% exposed")

    time = 0
    tmax = 1000
    dt = 0.1

    while time < tmax:
        # Use RK4 method to determine the new vector_population.
        k1 = derivs(time, vector_population)
        k2 = derivs(time + dt / 2, add(vector_population, mult(k1, 0.5)))
        k3 = derivs(time + dt / 2, add(vector_population, mult(k2, 0.5)))
        k4 = derivs(time + dt, add(vector_population, k3))
        vector_population = add(vector_population,
                                mult(add(mult(k1, 1 / 6), add(mult(k2, 1 / 3), add(mult(k3, 1 / 3), mult(k4, 1 / 6)))),
                                     dt))
        print(vector_population)

        # plot using vpython
        vpython.rate(50)  # fps of graph update
        s_plot.plot(time, vector_population[0])
        i_plot.plot(time, vector_population[1])
        r_plot.plot(time, vector_population[2])
        e_plot.plot(time, vector_population[3])

        time += dt


if __name__ == "__main__":
    main()
