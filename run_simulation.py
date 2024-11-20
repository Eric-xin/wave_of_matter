from calc import run_simulation

if __name__ == '__main__':
    # parameters: n, box_size, steps, initial_speed, particle_radius
    run_simulation(100, 100, 1000, 1000, 2)
    # the dpi of animation is set to 100 by default. it can be changed in the run_simulation function

    """
    Note that the you can also use the ParticleSimulation class directly to run the simulation:
    from calc import ParticleSimulation
    sim = ParticleSimulation(n=100, box_size=10, initial_speed=1, particle_radius=0.5, steps=1000)
    sim.run_simulation()
    """