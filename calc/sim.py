import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

# Constants
h = 6.62607015e-34  # Planck's constant (J·s)
c = 3.0e8           # Speed of light (m/s)
k_B = 1.380649e-23  # Boltzmann's constant (J/K)

class Particle:
    def __init__(self, position, velocity, mass=1.0, radius=0.5):
        self.position = np.array(position, dtype=float)
        self.velocity = np.array(velocity, dtype=float)
        self.mass = mass
        self.radius = radius
        self.update_kinetic_energy_and_temperature()

    def update_kinetic_energy_and_temperature(self):
        self.kinetic_energy = 0.5 * self.mass * np.linalg.norm(self.velocity)**2
        self.temperature = self.kinetic_energy / (1.5 * 1.38e-23)  # Boltzmann constant

def initialize_particles(n, box_size, initial_speed, particle_radius):
    particles = []
    for _ in range(n):
        position = np.random.rand(2) * box_size
        angle = np.random.rand() * 2 * np.pi
        velocity = initial_speed * np.array([np.cos(angle), np.sin(angle)])
        particles.append(Particle(position, velocity, radius=particle_radius))
    return particles

def update_positions(particles, dt, box_size):
    for particle in particles:
        particle.position += particle.velocity * dt
        # Handle wall collisions
        for i in range(2):
            if particle.position[i] - particle.radius <= 0 or particle.position[i] + particle.radius >= box_size:
                particle.velocity[i] *= -1
                particle.update_kinetic_energy_and_temperature()

def handle_collisions(particles):
    collision_count = 0
    for i, p1 in enumerate(particles):
        for j, p2 in enumerate(particles[i+1:], i+1):
            if np.linalg.norm(p1.position - p2.position) < p1.radius + p2.radius:
                # Inelastic collision: exchange velocities
                v1, v2 = p1.velocity, p2.velocity
                p1.velocity, p2.velocity = v2, v1
                
                # Introduce random fluctuations to simulate thermodynamic movement
                fluctuation = np.random.normal(0, 0.1, size=2)
                p1.velocity += fluctuation
                p2.velocity -= fluctuation
                
                # Update kinetic energy and temperature
                p1.update_kinetic_energy_and_temperature()
                p2.update_kinetic_energy_and_temperature()
                
                collision_count += 1
    return collision_count

def calculate_properties(particles):
    avg_speed = np.mean([np.linalg.norm(p.velocity) for p in particles])
    avg_temperature = np.mean([p.temperature for p in particles])
    entropy = np.sum([p.temperature * np.log(p.temperature) for p in particles])
    return avg_speed, avg_temperature, entropy

def de_broglie_wavelength(momentum):
    return h / momentum

def planck_law(wavelength, temperature):
    return (2 * h * c**2) / (wavelength**5 * (np.exp((h * c) / (wavelength * k_B * temperature)) - 1))

def animate(i, particles, scat, box_size, avg_speed_data, avg_temp_data, entropy_data, collision_freq_data, ax2):
    update_positions(particles, 0.01, box_size)
    collision_count = handle_collisions(particles)
    scat.set_offsets([p.position for p in particles])
    
    avg_speed, avg_temp, entropy = calculate_properties(particles)
    avg_speed_data.append(avg_speed)
    avg_temp_data.append(avg_temp)
    entropy_data.append(entropy)
    collision_freq_data.append(collision_count)
    
    momenta = [p.mass * np.linalg.norm(p.velocity) for p in particles]
    wavelengths = de_broglie_wavelength(np.array(momenta))
    
    ax2.cla()
    ax2.hist(wavelengths, bins=50, log=True)
    ax2.set_title("De Broglie Wavelengths")
    ax2.set_xlabel("Wavelength (m)")
    ax2.set_ylabel("Frequency (log scale)")
    ax2.grid(True)
    
    return scat,

def run_simulation(n, box_size, steps, initial_speed, particle_radius):
    particles = initialize_particles(n, box_size, initial_speed, particle_radius)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    ax1.set_xlim(0, box_size)
    ax1.set_ylim(0, box_size)
    max_length = box_size
    scat = ax1.scatter([p.position[0] for p in particles], [p.position[1] for p in particles], s=(particle_radius / max_length * 1000)**2)
    ax1.set_title("Particle Movement")
    ax1.set_xlabel("X Position")
    ax1.set_ylabel("Y Position")
    ax1.grid(True)

    # Plot Planck's black body radiation
    wavelengths_bb = np.linspace(1e-9, 3e-6, 500)  # Wavelength range (m)
    temperature = 5800  # Temperature (K)
    intensity = planck_law(wavelengths_bb, temperature)
    ax2.plot(wavelengths_bb * 1e9, intensity)  # Convert wavelength to nm for better visualization
    ax2.set_title("Planck's Black Body Radiation")
    ax2.set_xlabel("Wavelength (nm)")
    ax2.set_ylabel("Intensity")
    ax2.grid(True)

    cur_dir = os.getcwd()
    cur_dir = os.path.join(cur_dir, 'outputs')

    avg_speed_data, avg_temp_data, entropy_data, collision_freq_data = [], [], [], []

    ani = animation.FuncAnimation(fig, animate, fargs=(particles, scat, box_size, avg_speed_data, avg_temp_data, entropy_data, collision_freq_data, ax2), frames=steps, interval=50, blit=False)
    ani.save(cur_dir + '/particles_simulation.mp4', writer=animation.FFMpegWriter(fps=20, extra_args=['-vcodec', 'libx264', '-pix_fmt', 'yuv420p', '-b:v', '5000k']), dpi=100)
    
    # Plot average speed, temperature, and entropy curves
    plt.figure()
    plt.subplot(3, 1, 1)
    plt.plot(avg_speed_data)
    plt.title('Average Speed')

    plt.subplot(3, 1, 2)
    plt.plot(avg_temp_data)
    plt.title('Average Temperature')

    plt.subplot(3, 1, 3)
    plt.plot(entropy_data)
    plt.title('Entropy')

    plt.tight_layout()
    plt.savefig(cur_dir + '/average_speed_temperature_entropy.png')

    # Plot collision frequency curve
    plt.figure()
    plt.plot(collision_freq_data)
    plt.title('Collision Frequency')
    plt.savefig(cur_dir + '/collision_frequency.png')