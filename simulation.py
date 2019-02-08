import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
import math_functions

mf = math_functions.math_functions()

# Constants
au = 149597871.e3       # 1 AU (Astronomical Unit) = 149 597 871 km
G = 6.674e-11           # 6.674*10^(−11) N*(m/kg)^2

# Earth
m_earth = 5.97237e24    # 5.97237*10^(24) kg
radius_earth = 6371e3   # 6 371 km
area_earth = radius_earth**2 * np.pi
velocity_earth = 29.78 *1000   # 29.78 km/s

# Sun
m_sun = 333000*m_earth
radius_sun = 109*radius_earth
area_sun = (radius_sun/radius_earth)**2 * radius_earth

# Asteriod
m_asteriod = m_earth / 1.e10
radius_asteriod = (m_asteriod/m_earth)**(1/3) * radius_earth
area_asteriod = (radius_asteriod/radius_asteriod)**2 * radius_earth

# Time
minute = 60         # 60 s
hour = 60*minute    # 60 min
day = 24*hour       # 24 h
year = 365*day      # 365 days

# Initial conditions
r0_sun = [0,0,0]
r0_earth = [au,0,0]
r0_asteriod = [-au/2,au,0.2*au]
v0_sun = [0,0,0]
v0_earth = [0,velocity_earth,0]
v0_asteriod = [0,-velocity_earth,0]

class Simulation():

    def __init__(self, method = 'euler-cromer', adaptive_time_step = False):

        # ---------------------- Define variables

        # Constants
        au = 149597871.e3       # 1 AU (Astronomical Unit) = 149 597 871 km
        G = 6.674e-11           # 6.674*10^(−11) N*(m/kg)^2

        # Earth
        m_earth = 5.97237e24    # 5.97237*10^(24) kg
        radius_earth = 6371e3   # 6 371 km
        area_earth = radius_earth**2 * np.pi
        velocity_earth = 29.78 *1000   # 29.78 km/s

        # Sun
        m_sun = 333000*m_earth
        radius_sun = 109*radius_earth
        area_sun = (radius_sun/radius_earth)**2 * radius_earth

        # Asteriod
        m_asteriod = m_earth / 1.e10
        radius_asteriod = (m_asteriod/m_earth)**(1/3) * radius_earth
        area_asteriod = (radius_asteriod/radius_asteriod)**2 * radius_earth

        # Time
        minute = 60         # 60 s
        hour = 60*minute    # 60 min
        day = 24*hour       # 24 h
        year = 365*day      # 365 dys

        # Initial conditions
        r0_sun = [0,0,0]
        r0_earth = [au,0,0]
        r0_asteriod = [-au/2,au,0.2*au]
        v0_sun = [0,0,0]
        v0_earth = [0,velocity_earth,0]
        v0_asteriod = [0,-velocity_earth,0]

        # set up
        self.method = method
        self.methods_dict = {
            'euler': self.euler_step,
            'euler-cromer': self.euler_cromer_step,
            'verlet': self.verlet_step,
            'runge-kutta': self.runge_kutta_step
        }
        self.method_function = self.methods_dict[method]
        self.adaptive_time_step = adaptive_time_step

        # standard time data
        self.time_step = day
        self.t = 0
        self.t_end = year


        self.nr_of_objects = 0
        # Data to input, stored in following variables using add_object(...)
        self.name_array = []
        self.m_array = []
        self.r_array = []
        self.v_array = []
        self.color_data = []
        self.area_data = []

        # Calculated data
        self.x_data = []
        self.y_data = []
        self.z_data = []
        self.t_data = [self.t]
        self.distance_array = []    # between earth and asteriod
        self.kinetic_energy_data = []
        self.potential_energy_data = []
        self.total_energy_data = []

    def add_object(self, name, m,r,v,c='k',A=None):
        self.nr_of_objects += 1
        if name == None: name = str(self.nr_of_objects)
        self.name_array.append(name)
        self.m_array.append(m)
        self.r_array.append(r)
        self.v_array.append(v)

        self.x_data.append([r[0]])
        self.y_data.append([r[1]])
        self.z_data.append([r[2]])

        self.color_data.append(c)
        self.area_data.append(A)

    def execute_simulation(self):
        self.potential_energy_data.append(self._get_current_potential_energy())
        self.kinetic_energy_data.append(self._get_current_kinetic_energy())
        self.total_energy_data.append(self._get_current_total_energy())
        while self.t < self.t_end:
            self.method_function()
            self.store_current_iteration()

    def _get_current_potential_energy(self):
        """Returns the current potential energy of the system."""
        G = 6.674e-11           # 6.674*10^(−11) N*(m/kg)^2
        potential_energy = 0
        if self.nr_of_objects <= 1: return 0
        for i in range(1,self.nr_of_objects):
            for j in range(0, i-1):
                r_ij = mf.norm([x-y for x,y in zip(self.r_array[j], self.r_array[i])])
                potential_energy += -G*self.m_array[i]*self.m_array[j]/r_ij
        return potential_energy

    def _get_current_kinetic_energy(self):
        """Returns the current kinetic energy of the system."""
        kinetic_energy = 0
        for i in range(self.nr_of_objects):
            kinetic_energy += self.m_array[i]*mf.norm(self.v_array[i])**2
        return kinetic_energy/2


    def _get_current_total_energy(self):
        """Returns the current total energy of the system."""
        return self._get_current_potential_energy()+self._get_current_kinetic_energy()

    def store_current_iteration(self):
        self.t_data.append(self.t)
        self.potential_energy_data.append(self._get_current_potential_energy())
        self.kinetic_energy_data.append(self._get_current_kinetic_energy())
        self.total_energy_data.append(self._get_current_total_energy())
        for i in range(self.nr_of_objects):
            r = self.r_array[i]
            self.x_data[i].append(r[0])
            self.y_data[i].append(r[1])
            self.z_data[i].append(r[2])

    def euler_step(self):
        self.t += self.time_step
        # iterate through the objects
        for i in range(self.nr_of_objects):
            r_array = np.copy(self.r_array)
            j_range = [j for j in range(self.nr_of_objects) if j != i]
            # calculate acceleration
            a_i = [0,0,0]
            for j in j_range:
                a_i = mf.vector_add([
                    a_i,
                    self._get_a_ij(self.m_array[j],r_array[j],r_array[i])
                ])
            # update r
            self.r_array[i] = mf.vector_add([
                    self.r_array[i],
                    [v_i_i*self.time_step for v_i_i in self.v_array[i]]
            ])
            # update v
            self.v_array[i] = mf.vector_add([
                    self.v_array[i],
                    [a_i_i*self.time_step for a_i_i in a_i]
            ])

    def euler_cromer_step(self):
        self.t += self.time_step
        # iterate through the objects
        r_array = np.copy(self.r_array)
        for i in range(self.nr_of_objects):
            j_range = [j for j in range(self.nr_of_objects) if j != i]
            # calculate acceleration
            a_i = [0,0,0]
            for j in j_range:
                a_i = mf.vector_add([
                    a_i,
                    self._get_a_ij(self.m_array[j],r_array[j],r_array[i])
                ])
            # update v
            self.v_array[i] = mf.vector_add([
            self.v_array[i],
            [a_i_i*self.time_step for a_i_i in a_i]
            ])
            # update r
            self.r_array[i] = mf.vector_add([
                    self.r_array[i],
                    [v_i_i*self.time_step for v_i_i in self.v_array[i]]
            ])

    def verlet_step(self):
        print(self.t)
        self.t += self.time_step
        # iterate through the objects
        r_array = np.copy(self.r_array)
        for i in range(self.nr_of_objects):
            j_range = [j for j in range(self.nr_of_objects) if j != i]
            # calculate acceleration
            a_1 = [0,0,0]
            for j in j_range:
                a_1 = mf.vector_add([
                    a_1,
                    self._get_a_ij(self.m_array[j],r_array[j],r_array[i])
                ])
            # update r
            self.r_array[i] = mf.vector_add([
                    self.r_array[i],
                    [v_i_i*self.time_step for v_i_i in self.v_array[i]],
                    [1/2*a_1_i*(self.time_step**2) for a_1_i in a_1]
            ])
            # calculate acceleration at new r
            a_2 = [0,0,0]
            for j in j_range:
                a_2 = mf.vector_add([
                    a_2,
                    self._get_a_ij(self.m_array[j],r_array[j],self.r_array[i])
                ])
            # update v
            self.v_array[i] = mf.vector_add([
                    self.v_array[i],
                    [1/2*a_1_i*self.time_step for a_1_i in a_1],
                    [1/2*a_2_i*self.time_step for a_2_i in a_2]
            ])
            # uncomment when finding distance between asteriod and earth
            # R = [(self.r_array[3][i]-self.r_array[9][i]) for i in range(3)]
            # self.distance_array(mf.norm(R))


    def runge_kutta_step(self):
        self.t += self.time_step
        # iterate through the objects
        for i in range(self.nr_of_objects):
            # j_range is index for all objects other than i
            j_range = [j for j in range(self.nr_of_objects) if j != i]

            r = [x for x in self.r_array[i]]    # r(t)
            v = [x for x in self.v_array[i]]    # v(t)

            # a_1,b_1,...,a_4,b_4
            a_1 = [0,0,0]                       # = a(r(t))Dt
            for j in j_range:
                a_1 = mf.vector_add([
                    a_1,
                    self._get_a_ij(self.m_array[j],self.r_array[j],r)
                ])
            a_1 = mf.vector_scale(a_1,self.time_step)
            b_1 = mf.vector_scale(v,self.time_step) # = v(t)Dt

            r_2 = [r_i+b_1i/2 for r_i,b_1i in zip(r,b_1)] # = r(t)+b_1/2
            a_2 = [0,0,0]                       # = a(r_2)Dt
            for j in j_range:
                a_2 = mf.vector_add([
                    a_2,
                    self._get_a_ij(self.m_array[j],self.r_array[j],r)
                ])
            a_2 = mf.vector_scale(a_2,self.time_step)
            b_2 = [v_i+a_1i/2 for v_i,a_1i in zip(v,a_1)] # = v(t+a_1/2)Dt
            b_2 = mf.vector_scale(b_2, self.time_step)

            r_3 = [r_i+b_2i/2 for r_i,b_2i in zip(r,b_2)]   # = r(t)+b_2/2
            a_3 = [0,0,0]
            for j in j_range:
                a_3 = mf.vector_add([
                    a_3,
                    self._get_a_ij(self.m_array[j],self.r_array[j],r)
                ])
            a_3 = mf.vector_scale(a_3,self.time_step)
            b_3 = [v_i+a_2i/2 for v_i,a_2i in zip(v,a_2)]
            b_3 = mf.vector_scale(b_3,self.time_step)

            r_4 = [r_i+b_3i for r_i,b_3i in zip(r,b_3)]
            a_4 = [0,0,0]
            for j in j_range:
                a_4 = mf.vector_add([
                    a_4,
                    self._get_a_ij(self.m_array[j],self.r_array[j],r)
                ])
            a_4 = mf.vector_scale(a_4,self.time_step)
            b_4 = [v_i+a_3i/2 for v_i,a_3i in zip(v,a_3)]
            b_4 = mf.vector_scale(b_4,self.time_step)

            # update v
            self.v_array[i] = mf.vector_add([
                    self.v_array[i],
                    mf.vector_scale(
                    mf.vector_add([
                        a_1,
                        mf.vector_scale(a_2, 2),
                        mf.vector_scale(a_3, 2),
                        a_4
                    ])
                    , 1/6)
            ])

            # update r
            self.r_array[i] = mf.vector_add([
                    self.r_array[i],
                    mf.vector_scale(
                    mf.vector_add([
                        b_1,
                        mf.vector_scale(b_2, 2),
                        mf.vector_scale(b_3, 2),
                        b_4
                    ]),
                    1/6)
            ])

    def _get_a_ij(self, m2, r1, r2):
        """Takes mass m  """
        G = 6.674e-11           # 6.674*10^(−11) N*(m/kg)^2
        R = [(r1[i]-r2[i]) for i in range(3)]
        k = G*m2/(mf.norm(R)**3)
        return [k*R_i for R_i in R]

    def generate_animation(self, show = False, simulation_duration = 10, filename = None):
        nr_of_objects = self.nr_of_objects
        x_data = [[x/au for x in array] for array in self.x_data]
        y_data = [[y/au for y in array] for array in self.y_data]
        z_data = [[z/au for z in array] for array in self.z_data]
        data = []
        for i in range(nr_of_objects):
            data.append([x_data[i], y_data[i], z_data[i]])
        for i in range(nr_of_objects):
            data.append([x_data[i], y_data[i], z_data[i]])
        data = np.array(data)
        t_data = [t/day for t in self.t_data]
        color_array = self.color_data
        fig = plt.figure(figsize=(10,5), dpi = 100)
        ax1 = p3.Axes3D(fig)

        ax1.set_xlim(-5, 5)
        ax1.set_ylim(-5, 5)
        ax1.set_zlim(-5, 5)
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.set_zlabel('z')
        lines = []

        # Regulate real-time duration of Simulation
        simulation_duration = simulation_duration // 2 + 1
        number_of_frames = simulation_duration*1000 // 30
        i_step = len(t_data) / number_of_frames

        for i, dat in enumerate(data):  # grey lines
            if i < nr_of_objects:
                l = ax1.plot(dat[0,0:1],dat[1,0:1],dat[2,0:1],
                        linewidth = 0.5,
                        color = [0.5,0.5,0.5])[0]
                lines.append(l)
            else:
                l = ax1.plot([dat[0,1]],[dat[1,1]],[dat[2,1]],
                        markersize = 10,
                        marker = '.',
                        c = color_array[i%nr_of_objects])[0]
                lines.append(l)

        def animate(i, data_lines, lines):
            end_index = int(i*i_step)
            for j, line, data in zip(range(len(lines)), lines, data_lines):
                if j < nr_of_objects:
                    line.set_data(data[0:2,:end_index])
                    line.set_3d_properties(data[2,:end_index])
                else:
                    line.set_data(data[0:2,end_index])
                    line.set_3d_properties(data[2,end_index])
            return lines

        anim = animation.FuncAnimation(fig, animate, fargs=(data,lines),
                    interval=30, blit=True, frames = number_of_frames)
        if filename:
            anim.save(filename, dpi=100, writer='imagemagick')
            print('File succesfully saved.')
        if show: plt.show()

    def generate_energy_animation(self, show = False, simulation_duration = 10):
        """Not yet implemented"""
        pass
        # idea from: http://www.roboticslab.ca/wp-content/uploads/2012/11/robotics_lab_animation_example.txt

    def add_solar_system(self):
        """Sets initial conditions for simulation object .
        - color
        - plot_size
        -
        """

        self.add_object(
            name = 'Sun',
            m = 1.9885e30,      # mass of object in kg
            r = [0,0,0],        # initial radius from center of solar system
            v = [0,0,0],        # initial velocity
            c = [0.9,0.9,0],            # color in plot
            A = 696392000**2*np.pi  # cross-section area in m^2, for plot
        )
        self.add_object(
            name = 'Mecury',
            m = 3.3011e23,
            r = [57909050*1000,0,0],    # perihelion, m
            v = [0,47362,0],           # average orbital velocity
            c = [0.625,0.625,0.625],    # gray
            A = 2440000**2*np.pi
        )
        self.add_object(
            name = 'Venus',
            m = 4.8675e24,
            r = [0,108208000*1000,0],      # perihelion, km
            v = [-35020,0,0],               # average orbital velocity
            c = [160/256,160/256, 80/256],  # yellow
            A = np.pi*2440**2
        )
        self.add_object(
            name = 'Earth',
            m = 5.97237e24,
            r = [-149598023*1000,0,0],
            v = [0,-29780,0],
            c = [0.125, 0.400, 0.850],
            A = 6371000**2*np.pi
        )
        # uncomment below to add moon orbiting Earth
        # self.add_object( # TODO: Se till att denna cirkulerar jorden
        #     name = 'Moon',
        #     m = 7.342e22,
        #     r = [-149598023*1000,384399*1000,0],
        #     v = [1022,-29780,0],
        #     c = [0.8, 0.8, 0.8],
        #     A = 6371000**2*np.pi
        # )
        self.add_object(
            name = 'Mars',
            m = 6.4171e23,
            r = [0,-227939200*1000,0],   # semi-major axis
            v = [24007,0,0],
            c = [170/266,85/256,75/256],     # light red
            A = 3389999**2*np.pi
        )
        self.add_object(
            name = 'Jupiter',
            m = 1.8982e27,
            r = [778570000*1000,0,0],
            v = [0,13070,0],
            c = [210/256,165/256,85/256],
            A = 69911000**2*np.pi
        )
        self.add_object(
            name = 'Saturn',
            m = 5.6834e26,
            r = [0,1433530000*1000,0],
            v = [-9680,0,0],
            c = [210/256,200/256,35/256],
            A = 58232000**2*np.pi
        )
        self.add_object(
            name = 'Uranus',
            m = 8.6810e25,
            r = [-2875.04e9,0,0],
            v = [0,-6800,0],
            c = [110/256,160/256,222/256],
            A = 25362000**2*np.pi
        )
        self.add_object(
            name = 'Neptune',
            m = 1.02413e26,
            r = [0,-4.50e12,0],
            v = [5430,0,0],
            c = [47/256,53/256,222/256],
            A = 24622000**2*np.pi
        )
