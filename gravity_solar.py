import vpython

e_graph = vpython.gcurve(color=vpython.color.blue)


def gforce(p1, p2):
    # Calculate the gravitational force exerted on p1 by p2.
    G = 1  # Change to 6.67e-11 to use real-world values.
    # Calculate distance vector between p1 and p2.
    r_vec = p1.pos - p2.pos
    # Calculate magnitude of distance vector.
    r_mag = vpython.mag(r_vec)
    # Calcualte unit vector of distance vector.
    r_hat = r_vec / r_mag
    # Calculate force magnitude.
    force_mag = G * p1.mass * p2.mass / r_mag ** 2
    # Calculate force vector.
    force_vec = -force_mag * r_hat

    return force_vec


def main():

    star = vpython.sphere(pos=vpython.vector(0, 0, 0), radius=0.2, color=vpython.color.yellow, mass=1.0 * 1000, momentum=vpython.vector(0, 0, 0), make_trail=True)
    planet_one = vpython.sphere(pos=vpython.vector(1, 0, 0), radius=0.05, color=vpython.color.blue, mass=1, momentum=vpython.vector(0, 30, 0), make_trail=True)
    planet_two = vpython.sphere(pos=vpython.vector(0, 1, 0), radius=0.1, color=vpython.color.red, mass=1, momentum=vpython.vector(- 35, 0, 0), make_trail=True)
    planet_three = vpython.sphere(pos=vpython.vector(0, 8, 0), radius=0.1, color=vpython.color.red, mass=100, momentum=vpython.vector(- 100, 0, 0), make_trail=True)


    dt = 0.0001
    t = 0
    while True:
        vpython.rate(500)

        # Calculate forces.
        star.force = gforce(star, planet_one) + gforce(star, planet_two) + gforce(star, planet_three)
        planet_one.force = gforce(planet_one, star) + gforce(planet_one, planet_two) + gforce(planet_one, planet_three)
        planet_two.force = gforce(planet_two, star) + gforce(planet_two, planet_one) + gforce(planet_two, planet_three)
        planet_three.force = gforce(planet_three, star) + gforce(planet_three, planet_one) + gforce(planet_three, planet_two)


        # Update momenta.
        star.momentum = star.momentum + star.force * dt
        planet_one.momentum = planet_one.momentum + planet_one.force * dt
        planet_two.momentum = planet_two.momentum + planet_two.force * dt
        planet_three.momentum = planet_three.momentum + planet_three.force * dt

        # Update positions.
        star.pos = star.pos + star.momentum / star.mass * dt
        planet_one.pos = planet_one.pos + planet_one.momentum / planet_one.mass * dt
        planet_two.pos = planet_two.pos + planet_two.momentum / planet_two.mass * dt
        planet_three.pos = planet_three.pos + planet_three.momentum / planet_two.mass * dt

        t = t + dt


if __name__ == "__main__":
    main()
