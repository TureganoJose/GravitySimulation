import numpy as np  # need arrays
import matplotlib.pyplot as plt  # need plots

# By default all units are SI
# Positive -> Upward (gravity always negative)

# Main constants
simulation_time = 10  # [s]
PI = 3.141592
G = 6.67408e-11  # Gravitational constant[m3 kg-1 s-2]
dt = 0.001  # step time[s]


class Fluid:
   def __init__(self, density):
       self.density = density


class Planet:
   def __init__(self, radius, mass, gravity_acc):
       self.radius = radius
       self.mass = mass
       self.g = gravity_acc


class Object:
   def __init__(self, mass, drag_coeff, radius, state, fluid, planet):
       if mass < 0.000001:
           print('Warning: very small or close to zero mass')
       self._mass = mass
       self._drag_coeff = drag_coeff
       self._radius = radius
       self._state = state  # 0) Position 1) Velocity
       # self.__derivative = derivative  # 0) Velocity 1) Acceleration
       self._fluid = fluid
       self._planet = planet

   def drag_force(self, speed):
       return - np.sign(speed) * self._drag_coeff * PI * self._radius * self._radius * self._fluid.density * speed * speed

   def gravity(self):
       return - self._planet.g * self._mass
   def gravity(self, height_pos):
       return -G * self._planet.mass * self._mass/((self._planet.radius+height_pos)*(self._planet.radius+height_pos))
   # WATCH OUT FOR DIVIDE BY ZERO

   def calculate_acceleration(self, input_state):
       # Equation of motion
       return np.array([input_state[1], (self.drag_force( input_state[1]) + self.gravity(input_state[0])) / self._mass])  # 0) Velocity 1) Acceleration

   def first_derivative(self):
       derivative_state = self.calculate_acceleration(self._state)
       return derivative_state

   def next_derivative(self, derivative, h):
       # temp_state = np.array([0.0, 0.0])  # temp state vector: 0) Position 1) Velocity
       temp_state = self._state + derivative * h
       return self.calculate_acceleration(temp_state)

   def integrator_rk4(self):
       k1 = self.first_derivative()
       k2 = self.next_derivative( k1, dt * 0.5)
       k3 = self.next_derivative( k2, dt * 0.5)
       k4 = self.next_derivative( k3, dt)
       self._state = self._state + dt* (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

   def integrator_euler(self):
       self._state = self._state + dt * self.calculate_acceleration(self._state)


# Inputs
drop_height = 100.0 # [m]
initial_velocity = 10.0 # [m s-1] Positive -> Upward
air = Fluid(density=1.21)
earth = Planet(radius=6371e3, mass=5.972e24, gravity_acc=9.81)
ball = Object(mass=2, drag_coeff=1.0, radius=0.1, state=np.array([drop_height, initial_velocity]), fluid=air, planet=earth)

# logging variables
position = []
velocity = []
acceleration =[]
t = []
position_check=[]
velocity_check=[]
acceleration_check=[]

# Main loop
time = 0
while time < simulation_time:
   # Update time step
   ball.integrator_rk4()

   # logging
   position.append(ball._state[0])
   velocity.append(ball._state[1])
   [ speed, acc] = ball.calculate_acceleration(ball._state)
   acceleration.append(acc)

   # logging check physics (no drag)
   position_check.append(drop_height - earth.g * time*time/2)
   velocity_check.append( -np.sqrt((ball._planet.g * ball._mass)/(ball._drag_coeff * PI * ball._radius * ball._radius * ball._fluid.density))) # terminal speed as reference
   acceleration_check.append(-earth.g)
   t.append(time)
   time = time + dt

# Plotting
f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
ax1.plot(t,position, 'r', label='Constant Gravity')
ax1.plot(t, position_check, 'k', label='No drag')
ax2.plot(t,velocity, 'r', label='Constant Gravity')
ax2.plot(t,velocity_check, 'k', label='Terminal velocity (Drag)')
ax3.plot(t,acceleration, 'r', label='Constant Gravity')
ax3.plot(t,acceleration_check, 'k', label='No drag')

ax1.set(ylabel='Position [m]')
ax2.set(ylabel='Speed [m/s]')
ax3.set(ylabel='Acceleration [m/s2]')
ax3.set(xlabel='time[s]')
ax1.legend(loc="upper right")
ax2.legend(loc="upper right")
ax3.legend(loc="upper right")
plt.show()
