# Simulation of the solar system

This project simulates the solar system using four different numerical methods.

## Results

The following are some example animations that can be generated using this libary.  This is the result when a **light** object (*black*) passes close to the earth (*blue*).

![light](light-object.gif)

This is the result when a **heavy** object (*black*) passes close to the earth (*blue*).

![heavy](heavy-object.gif)

## Using the *simulation.py* module

### Installation

Please place the `simulation.py` and `math_functions.py` in a repository where you would like to use the module. You also need the following Python libraries:

- `matplotlib`
- `numpy`

### Usefull functions

The simulation module solves the equations of [Newton's law of gravitation](https://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion#Newton's_law_of_gravitation) with stated initial conditions. With the library a simulation object is initiated with:

```python
# methods: 'euler', 'euler-cromer', 'verlet' or 'runge-kutta'

simulation = Simulation('euler')
simulation.time_step = 10*day # in seconds, standard: 1 hour
simulation.time_end = 100*year  # in seconds, standard: 1 year
```

The function  `simulation.add_solar_system()` adds all planets in the solar system. With the `simulation.add_object()` function objects can be added to the simulation.  For example adding an asteriod with collision course to Earth can be done with

```python
simulation.add_object(
    name = 'Asteriod',	# name of object
    m = 9.393e20,		# mass of object, kg
    r = [414010000*1000,0,0],	# initial displacement vector, m
    v = [-17905/4*1.08497526,17905*0.5,1000],	# initial velocity vector, m/s
    c = [50/256,60/256,60/256],	# color of object for plot, optional
    A = 473000**2*np.pi		# surface area of object for plot, optional 
)
```

The simulation is then executed and a and a 3D-rendering of simulation can be generated with the script

```python
simulation.execute_simulation()
simulation.generate_animation(
	show=True,	# shows animation in new window
	filename='simulation_movie.gif',	# optional
	simulation_duration = 10 		# approximate duration for the simulation, s, default is 10
)
```

## Example 

In the `main.py` the following function adds an asteriod to the simulation:

```python
def add_asteriod(simulation, m=9.393e20, vz=0):
    """Adds an asteriod to a simulation object this mass m and vertical velocity vz.
    """
    simulation.add_object(
        name = 'Asteriod',
        m = m,
        r = [414010000*1000,0,0],
        # v = [-17905/4*1.08497526,17905*0.5,vz],
        v = [-5000,5000,vz],
        c = [50/256,60/256,60/256],
        A = 473000**2*np.pi
    )
```

Using this function an animation can be generated using the following function:

```python
def generate_animation():
    # uncomment below to change mass of asteroid
    # m = None # ceres
    # m = m_earth
     m = 317.8*m_earth # jupiter
    # m = m_sun

    simulation = Simulation('verlet')
    simulation.add_solar_system()
    add_asteriod(simulation, m)
    simulation.time_step = day
    simulation.t_end = year*2
    simulation.execute_simulation()
    simulation.generate_animation(simulation_duration = 10, show=True)
```
