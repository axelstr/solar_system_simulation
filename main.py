# Asteroid tragectory simulation
# SI1336
# Axel Stromberg
# axelstr@kth.se
#
# Project:
# Simulate our solar system, add an asteroid and make it pass very close
# (such that it’s direction changes signiﬁcantly) by a planet, you can also vary
# the mass of the asteroid from realistic to planet like.
#
# SI-units:
# - length: m
# - time:   s
# - mass:   kg



# ---------------------- Import libraries
# To save animation
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.pylab import *
# from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
# import math as m
import datetime
date = str(datetime.datetime.now()).replace(':','.')[2:19]
dated_filename = 'assets/archive/'+date

from simulation import Simulation
import math_functions as mf

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['font.size'] = 12


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

# ---------------------- Initial conditions

def set_initial_conditions_only_earth(simulation):
    """Sets initial conditions for simulation object .
    - color
    - plot_size
    -
    """
    simulation.time_step = day
    simulation.t_end = year

    simulation.add_object(
        name = 'Sun',
        m = 1.9885e30,      # mass of object in kg
        r = [0,0,0],        # initial radius from center of solar system
        v = [0,0,0],        # initial velocity
        c = [0.9,0.9,0],            # color in plot
        A = 696392000**2*np.pi  # cross-section area in m^2, for plot
    )
    simulation.add_object(
        name = 'Earth',
        m = 5.97237e24,
        r = [-147095000*1000,0,0],
        v = [0,-29780,0],
        c = [0.125, 0.400, 0.850],
        A = 6371000**2*np.pi
    )

def add_asteriod(simulation, m=9.393e20, vz=0):
    """Sets initial conditions for simulation object .
    - color
    - plot_size
    -
    """
    # simulation.add_object(
    #     name = 'Asteriod',
    #     m = m,
    #     r = [414010000*1000,0,0],
    #     # v = [-17905/4*1.08497526,17905*0.5,vz],
    #     v = [-17905/4*1.080055,17905*0.48640,vz],
    #     c = [50/256,60/256,60/256],
    #     A = 473000**2*np.pi
    # )
    simulation.add_object(
        name = 'Asteriod',
        m = m,
        r = [414010000*1000,0,0],
        # v = [-17905/4*1.08497526,17905*0.5,vz],
        v = [-5000,5000,vz],
        c = [50/256,60/256,60/256],
        A = 473000**2*np.pi
    )

def add_asteriod_and_sun(simulation,m=9.393e20):
    simulation.add_object(
        name = 'Sun',
        m = 1.9885e30,      # mass of object in kg
        r = [0,1,0],        # initial radius from center of solar system
        v = [0,1,0],        # initial velocity
        c = [0.9,0.9,0],            # color in plot
        A = 696392000**2*np.pi  # cross-section area in m^2, for plot
    )
    simulation.add_object(
        name = 'Asteriod',
        m = m,
        r = [414010000*1000,100,0],
        # v = [-17905/4*1.08497526,17905*0.5,vz],
        v = [0,10000,0],
        c = [50/256,60/256,60/256],
        A = 473000**2*np.pi
    )

# ---------------------- Main calculation

def get_initial_positions_plot():
    simulation = Simulation('euler-cromer')
    simulation.add_solar_system()
    x_data = [[x/au for x in array] for array in simulation.x_data]
    y_data = [[y/au for y in array] for array in simulation.y_data]
    color_array = simulation.color_data
    for i in [x for x in range(simulation.nr_of_objects) if x!=4]:
         plt.scatter(x_data[i],y_data[i],
                s = 50,
                marker = '.',
                c = [color_array[i%simulation.nr_of_objects]],
                zorder = 2)
    # plt.axis([-10,10,-10,10])
    plt.axis('equal')
    plt.grid(True, linestyle = ':', markevery=10, zorder = 1)
    plt.legend([name for name in simulation.name_array if name!='Moon'], loc='lower left')
    plt.title('Initial positions of the planets\nrotating counter-clockwise '+r'$\circlearrowleft$')
    plt.xlabel('x/[au]')
    plt.ylabel('y/[au]')
    date = str(datetime.datetime.now()).replace(':','.')[2:19]
    plt.savefig('assets/archive/plot %s.png'%(date), dpi = 300)
    plt.show()

def generate_animation():
    # uncomment below to change mass of asteroid
    # m = None # ceres
    # m = m_earth
    # m = 317.8*m_earth # jupiter
    # m = m_sun

    simulation = Simulation('verlet')
    simulation.add_solar_system()
    add_asteriod(simulation, m)
    # add_asteriod_and_sun(simulation)
    simulation.time_step = day
    simulation.t_end = year*2
    simulation.execute_simulation()
    simulation.generate_animation(simulation_duration = 10, show=True) #,filename=dated_filename+'.gif')
    # print(min(simulation.distance_array)/384399*1000)

def energy_plot():
    methods = ['euler', 'euler-cromer', 'verlet', 'runge-kutta']
    titles = ['Euler', 'Euler-Cromer','Verlet','Runge-Kutta']
    common_linewidth = 1
    c1,c2,c3 = [40/256,40/256,251/256], [59/256,122/256,87/256], [204/256,0,0]
    for i in range(4):
        method = methods[i]

        # simulation
        simulation = Simulation(method)
        simulation.add_solar_system()
        simulation.time_step = day*10
        # simulation.time_step = day*1
        simulation.t_end = 100*year
        simulation.execute_simulation()
        print('Currently at ' + method)

        # plot
        initial_total_energy = np.abs(simulation.total_energy_data[0])
        plt.subplot(2,2,i+1)
        plt.plot([t/year for t in simulation.t_data], [E/initial_total_energy for E in simulation.total_energy_data],
            color = c1,
            linewidth = common_linewidth,
            zorder = 2)
        plt.plot([t/year for t in simulation.t_data], [E/initial_total_energy for E in simulation.potential_energy_data],
            color = c2,
            linewidth = common_linewidth)
        plt.plot([t/year for t in simulation.t_data], [E/initial_total_energy for E in simulation.kinetic_energy_data],
            color = c3,
            linewidth = common_linewidth)
        plt.plot([0, simulation.t_data[-1]/year], [-1,-1],
            color = 'k',
            linewidth = common_linewidth,
            linestyle = ':',
            zorder = 2)
        plt.title(titles[i])
        if i==2:
            plt.ylabel('Energy / [initial energy]')
            plt.xlabel('time / [years]')
        if i==3:
            plt.legend([r'$E_{tot}$',r'$E_p$',r'$E_k$'])
    plt.tight_layout()
    plt.savefig('assets/archive/plot %s.png'%(date), dpi = 300)
    plt.show()

if __name__ == '__main__':
    # get_initial_positions_plot()
    # energy_plot()
    generate_animation()




# euler_simulation.generate_animation(
# 	show=True,	# shows animation in new window
#     # filename='simulation_movie.mp4',	# optional
# 	simulation_duration = 5 		# approximate duration for the simulation, s
# )
# print(min(euler_simulation.TEMParray))

#---------------------- Analysis

# ---------------------- Plot

# x_earth = euler_simulation.x_data[1]
# y_earth = euler_simulation.y_data[1]
# x_asteriod = euler_simulation.x_data[2]
# y_asteriod = euler_simulation.y_data[2]

# print('x_earth',x_earth)
# print('y_earth',y_earth)
# plt.plot(x_earth, y_earth)
# plt.plot(x_asteriod, y_asteriod)
# plt.show()

# -------------------- Animate





# # line1_2 = ax1.plot(0,0)
#
# def update_lines(i):
#     line1_1.set_data(data[:i])
#     return line1_1,
#
# fig, [ax1, ax2] = plt.subplots(1,2)
# data = [x_earth, y_earth]
# line1_1 = ax1.plot([],[],'r-')
#
# ani = animation.FuncAnimation(
#         fig,
#         update_lines,
#         interval = 2,
#         blit = True,
#         fargs = (data,line1_1))
#
# plt.show()


# def update_line(num, data, line):
#     line.set_data(data[...,:num])
#     return line,
#
# fig1 = plt.figure()
#
# # data = [x_earth, y_earth]
# data = np.random.rand(2, 25)
# l, = plt.plot([], [], 'r-')
# # plt.xlim(-2*au, 2*au)
# # plt.ylim(-2*au, 2*au)
# plt.xlabel('x')
# plt.title('test')
# line_ani = animation.FuncAnimation(fig1, update_line, 25, fargs=(data, l),
#     interval=50, blit=True)
# plt.show()
#line_ani.save('lines.mp4')'


# # WORKING:
# def update_line(num, data, line):
#     line.set_data(data[...,:num])
#     return line,
#
# fig1 = plt.figure()
#
# data = np.array([x_earth, y_earth])
# l, = plt.plot([], [], 'r-')
# plt.xlim(-2*au, 2*au)
# plt.ylim(-2*au, 2*au)
# plt.xlabel('x')
# plt.title('test')
# line_ani = animation.FuncAnimation(fig1, update_line, 25, fargs=(data, l),
#     interval=50, blit=True)
# plt.show()
# #line_ani.save('lines.mp4')
