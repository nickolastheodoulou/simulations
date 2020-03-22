import vpython

# GlowScript 2.9 VPython
# https://www.idmod.org/docs/hiv/model-seir.html
# https://en.wikipedia.org/wiki/Semi-implicit_Euler_method
# https://sci-hub.tw/10.1007/s10483-007-0914-x


def derivs(t, population, parameters):
    # parameters
    b = parameters['birth_rate']  # birth rate
    mu = parameters['death_rate']  # death rate
    beta = parameters['rate_of_transmission']  # rate of transmission
    alpha = parameters['rate_of_incubation']  # the rate of latent individuals becoming infectious (average duration of incubation is 1/\sigma)
    gamma = parameters['rate_of_recovery']  # the rate of cure of disease
    sigma = parameters['rate_of_loss_of_immunity']  # the rate at which recovered individuals return to the susceptible statue due to loss of immunity.
    nu = parameters['rate_of_vaccination']  # the ratio of the number of individuals who received the vaccine
    q = parameters['rate_of_treatment']  # giving treatment

    s = population[0]  # number susceptible: Not yet infected.
    e = population[1]  # number exposed
    i = population[2]  # number infected: Ongoing disease.
    r = population[3]  # number recovered: Gained immunity, isolated, or dead.
    total_population = s + e + i + r

    ds_dt = (b * total_population) + (sigma * r) - (beta * s * i / total_population) - (nu * s) - (mu * s)
    de_dt = (beta * s * i / total_population) - (alpha * e) - (mu * e)
    di_dt = (alpha * e) + (q * i) - (gamma * i) - (mu * i)
    dr_dt = (gamma * i) + (nu * s) - (sigma * r) - (mu * r)

    return [ds_dt, de_dt, di_dt, dr_dt]


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

    parameters = {
        "birth_rate": 0.00,
        "death_rate": 0.00,
        "rate_of_transmission": 0.96,
        "rate_of_incubation": 1,
        "rate_of_recovery": 0.90,
        "rate_of_loss_of_immunity": 0.003,
        "rate_of_vaccination": 0,
        "rate_of_treatment": 0

    }

    pop_susceptible = 50000000  # number susceptible: Not yet infected.
    pop_infectious = 500  # number infected: Ongoing disease.
    pop_recovered = 100  # number recovered: Gained immunity, isolated, or dead.
    pop_exposed = 5000000  # number exposed
    vector_population = [pop_susceptible, pop_exposed, pop_infectious, pop_recovered]

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
        k1 = derivs(time, vector_population, parameters)
        k2 = derivs(time + dt / 2, add(vector_population, mult(k1, 0.5)), parameters)
        k3 = derivs(time + dt / 2, add(vector_population, mult(k2, 0.5)), parameters)
        k4 = derivs(time + dt, add(vector_population, k3), parameters)
        vector_population = add(vector_population, mult(add(mult(k1, 1 / 6), add(mult(k2, 1 / 3), add(mult(k3, 1 / 3), mult(k4, 1 / 6)))), dt))
        print(vector_population)

        # plot using vpython
        vpython.rate(50)  # fps of graph up§te
        s_plot.plot(time, vector_population[0])
        e_plot.plot(time, vector_population[1])
        i_plot.plot(time, vector_population[2])
        r_plot.plot(time, vector_population[3])


        time += dt


if __name__ == "__main__":
    main()
