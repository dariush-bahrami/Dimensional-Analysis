```python
import sympy
sympy.init_printing()

from dimensional_analysis import (parameter, dimensional_analysis,
                                  solve_from_dimensional_analysis,
                                  standard_dimensional_analysis,
                                  solve_from_standard_dimensional_analysis)
```

# Manual Dimensional Analysis


```python
density = parameter('density', 'kg m^-3', '\\rho')
viscosity = parameter('viscosity', 'kg m^-1 s^-1', '\mu')
velocity = parameter('velocity', 'm s^-1', 'u')
diameter = parameter('diameter', 'm', 'D')
```


```python
dimensional_analysis(density, viscosity, velocity, diameter)
```



$$
\displaystyle \left\{ \Pi_{0} : \frac{\rho D u}{\mu}\right\}
$$



```python
solve_from_dimensional_analysis(density, viscosity, velocity, diameter, target_parameter=velocity)
```



$$
\displaystyle \left\{ \Pi_{0} : \left[ \frac{\Pi_{0} \mu}{\rho D}\right]\right\}
$$


# Automated Standard Dimensional Analysis

## Kinetic Energy


```python
standard_dimensional_analysis('kinetic_energy mass velocity')
```



$$
\displaystyle \left\{ \Pi_{0} : \frac{\sqrt{m} v}{\sqrt{KE}}\right\}
$$



```python
solve_from_standard_dimensional_analysis('kinetic_energy mass velocity', 'mass')
```



$$
\displaystyle \left\{ \Pi_{0} : \left[ \frac{\Pi_{0}^{2} KE}{v^{2}}\right]\right\}
$$


## Fluid Flow


```python
standard_dimensional_analysis('density viscosity velocity diameter')
```



$$
\displaystyle \left\{ \Pi_{0} : \frac{\rho D v}{\mu}\right\}
$$


## Pendulum


```python
standard_dimensional_analysis('period g length')
```



$$
\displaystyle \left\{ \Pi_{0} : \frac{L}{T^{2} g}\right\}
$$


## Fluid Static


```python
standard_dimensional_analysis('pressure density g height')
```



$$
\displaystyle \left\{ \Pi_{0} : \frac{\rho h g}{P}\right\}
$$


## Wave Equation


```python
standard_dimensional_analysis('wave_length velocity period')
```



$$
\displaystyle \left\{ \Pi_{0} : \frac{T v}{\lambda}\right\}
$$


## Magnetic Force


```python
standard_dimensional_analysis('magnetic_field force electric_current length')
```



$$
\displaystyle \left\{ \Pi_{0} : \frac{I L B}{F}\right\}
$$



```python
solve_from_standard_dimensional_analysis('magnetic_field force electric_current length', 'force')
```



$$
\displaystyle \left\{ \Pi_{0} : \left[ \frac{I L B}{\Pi_{0}}\right]\right\}
$$
