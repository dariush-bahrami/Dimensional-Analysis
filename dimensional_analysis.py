def raise_unit_2_power(unit_string, power):
    """Raise unit string to power

    Parameters
    ----------
    unit_string : str
    power : float

    Returns
    -------
    str
    """
    unit_list = unit_string.split()
    new_unit_list = []
    if power == 1:
        return unit_string
    for unit in unit_list:
        if '^' in unit:
            unit_part = unit.split('^')[0]
            power_part = float(unit.split('^')[1])
            power_part *= power

            if power_part - int(power_part) == 0:
                power_part = int(power_part)

            new_unit_list.append('^'.join([unit_part, str(power_part)]))
        else:
            unit_part = unit
            power_part = power
            if power_part - int(power_part) == 0:
                power_part = int(power_part)
            new_unit_list.append('^'.join([unit_part, str(power_part)]))

    return ' '.join(new_unit_list)


def si_derived_unit_equivalent(unit_string: str) -> str:
    """Convert SI derived units to fundamental
    for example: 'N' (newton) wil became 'kg m s^-2'

    Parameters
    ----------
    unit_string : str
        
    Returns
    -------
    str        
    """
    assert ' ' not in unit_string, 'Invalid space character'

    si_derived_units = {
        'N': 'kg m s^-2',
        'J': 'kg m^2 s^-2',
        'C': 'A s',
        'T': 'kg s^-2 A^-1',
        'Pa': 'kg m^-1 s^-2',
        'W': 'kg m^2 s^-3',
        'Hz': 's^-1',
        'V': 'kg m^2 s^-3 A^-1',
        'F': 'kg^-1 m^-2 s^4 A^2',
        'Wb': 'kg m^2 s^-2 A^-1',
        'H': 'kg m^2 s^-2 A^-2',
        'ohm': 'kg m^2 s^-3 A^-2',
        'rad': '',
        'sr': ''
    }

    if '^' in unit_string:
        unit_part = unit_string.split('^')[0]
        power_part = float(unit_string.split('^')[1])
    else:
        unit_part = unit_string
        power_part = 1

    if unit_part in si_derived_units:
        unit_part = si_derived_units[unit_part]
        return raise_unit_2_power(unit_part, power_part)

    else:
        return unit_string


def simplify_derived_units(unit_string: str) -> str:
    """Batch conversion of SI derived units

    Parameters
    ----------
    unit_string : str
        
    Returns
    -------
    str
    """
    unit_list = []
    for unit in unit_string.split(' '):
        unit_list.append(si_derived_unit_equivalent(unit))
    return ' '.join(unit_list)


def si_parser(unit_string: str) -> dict:
    """Convert SI units to fundamental system-independent dimensions.
    for example: 'm' will became sympy.physics.units.length

    Parameters
    ----------
    unit_string : str
        
    Returns
    -------
    dict
        dictionary of basic dimensions with coressponding power        
    """
    from sympy.physics import units

    si_system = {
        'm': units.length,
        's': units.time,
        'k': units.temperature,
        'kg': units.mass,
        'mol': units.amount_of_substance,
        'A': units.current,
        'cd': units.luminous_intensity
    }
    unit_string = simplify_derived_units(unit_string)

    result_unit = 1
    for unit in unit_string.split(' '):

        if '^' in unit:
            unit_part = unit.split('^')[0]
            power_part = int(unit.split('^')[1])
        else:
            unit_part = unit
            power_part = 1

        result_unit *= (si_system[unit_part]**power_part)

    return result_unit


def parameter(name: str, unit: str, latex_repr: str):
    """Simple function for creating sympy Quantity objects

    Parameters
    ----------
    name : str
        name of the physical quantity
    unit : str
        unit of quantity; for example unit of velocity is 'm s^-1'
    latex_repr : str
        latext representation of quantity; for example for density this
        parameter is: '\\rho' and for viscosity is '\mu'

    Returns
    -------
    sympy Quantity object
    """
    from sympy.physics.units import Quantity
    from sympy.physics.units.systems import SI

    parameter = Quantity(name, latex_repr=latex_repr)
    SI.set_quantity_dimension(parameter, si_parser(unit))
    return parameter


class DimensionalAnalysis:
    """Python class for Buckingham Theorem 
    """
    def __init__(self, parameters: list):
        """
        Parameters
        ----------
        parameters : list
            list of sympy Quantity objects
        """
        self.parameters = parameters

    @property
    def fundamental_dimensions(self):
        from sympy.physics.units.systems.si import dimsys_SI
        dimensions = set()
        for parameter in self.parameters:
            dimension_dict = dimsys_SI.get_dimensional_dependencies(
                parameter.dimension)
            for dimension in dimension_dict:
                dimensions.add(dimension)
        return dimensions

    @property
    def dimension_matrix(self):
        import sympy
        from sympy.physics.units.systems.si import dimsys_SI
        matrix = sympy.zeros(len(self.fundamental_dimensions),
                             len(self.parameters))
        for i, dimension in enumerate(self.fundamental_dimensions):
            for j, parameter in enumerate(self.parameters):
                dimension_dict = dimsys_SI.get_dimensional_dependencies(
                    parameter.dimension)
                if dimension in dimension_dict:
                    matrix[i, j] = dimension_dict[dimension]
                else:
                    matrix[i, j] = 0

        return matrix

    @property
    def dimensionless_parameters(self):
        import sympy
        dimesionless_dict = dict()
        nullspace = self.dimension_matrix.nullspace()
        for i, vector in enumerate(nullspace):
            d = 1
            for j, power in enumerate(vector):
                if power - int(power) == 0:
                    power = int(
                        power)  # Prefer using integer powers if possible
                d *= (self.parameters[j]**power)
                dimesionless_dict[sympy.symbols(f'Pi_{i}')] = d

        return dimesionless_dict

    def solve_for(self, parameter) -> dict:
        """Solve result from dimensional analysis for selected parameter

        Parameters
        ----------
        parameter : sympy Quantity object
            selected parameter for solving result from dimensional analysis

        Returns
        -------
        dict
            dictionary of solutions for each dimensionless parameter
        """
        import sympy
        solution_dict = dict()
        dimensionless_dict = self.dimensionless_parameters
        for d in dimensionless_dict:
            solution = sympy.solve(sympy.Eq(dimensionless_dict[d], d),
                                   parameter)
            solution_dict[d] = solution
        return solution_dict


def dimensional_analysis(*args) -> dict:
    """Performing dimensional analysis on input Quantiy objects

    Returns
    -------
    dict
        dictionary of dimensionless parameters
    """
    analyze = DimensionalAnalysis([*args])
    return analyze.dimensionless_parameters


def solve_from_dimensional_analysis(*args, target_parameter) -> dict:
    """Performing dimensionall analysis and solving results for specific
    parameter

    Parameters
    ----------
    target_parameter : sympy Quantity object
        selected parameter for solving result from dimensional analysis

    Returns
    -------
    dict
        dictionary of solutions for each dimensionless parameter
    """
    analyze = DimensionalAnalysis([*args])
    return analyze.solve_for(target_parameter)


def standard_parameters(parameters_string: str) -> list:
    """Collection of standard physical parameters

    Parameters
    ----------
    parameters_string : str
        name of parameters which will be returned as sympy Quantity object; for
        example: 'density viscosity length'

    Returns
    -------
    list
        list of sympy Quantity objects
    """
    parameter_list = parameters_string.split(' ')

    # 1D Length related parameters
    length = parameter('Length', 'm', 'L')
    width = parameter('Width', 'm', 'L')
    height = parameter('Height', 'm', 'h')
    diameter = parameter('Diameter', 'm', 'D')
    radius = parameter('Radius', 'm', 'r')
    amplitude = parameter('Amplitude', 'm', 'A')
    thickness = parameter('Thickness', 'm', 't')

    # 2D Length related parameters
    area = parameter('Area', 'm^2', 'A')

    # 3D Length related parameters
    volume = parameter('Volume', 'm^2', 'V')

    # Mechanic related parameters
    mass = parameter('Mass', 'kg', 'm')
    velocity = parameter('Velocity', 'm s^-1', 'v')
    speed = parameter('Speed', 'm s^-1', 'v')
    angular_velocity = parameter('Angular Velocity', 's^-1', '\omega')
    acceleration = parameter('Acceleration', 'm s^-2', 'a')
    g = parameter('g', 'm s^-2', 'g')
    force = parameter('Force', 'N', 'F')
    momentum = parameter('Momentum', 'kg m s^-1', 'p')
    period = parameter('Period', 's', 'T')
    frequency = parameter('Frequency', 's^-1', 'f')
    energy = parameter('Energy', 'J', 'E')
    work = parameter('Work', 'J', 'W')
    potential_energy = parameter('Potential Energy', 'J', 'PE')
    kinetic_energy = parameter('Kinetic Energy', 'J', 'KE')
    tension = parameter('Tension', 'N', 's')
    linear_density = parameter('Linear Density', 'kg m^-1', '\\rho')
    stress = parameter('Stress', 'N m^-2', '\sigma')
    power = parameter('Power', 'J s^-1', 'P')

    # Fluid mechanic parameters
    density = parameter('Density', 'kg m^-3', '\\rho')
    viscosity = parameter('Viscosity', 'kg m^-1 s^-1', '\mu')
    pressure = parameter('Pressure', 'Pa', 'P')
    temperature = parameter('Temperature', 'k', 'T')
    heat = parameter('Heat', 'J', 'Q')

    # Wave related parameters
    wave_length = parameter('Wave Length', 'm', '\lambda')

    # Electric related parameters
    electric_current = parameter('Electric Current', 'A', 'I')
    electric_voltage = parameter('Electric Voltage', 'V', 'V')
    electric_resistance = parameter('Electric Resistance', 'ohm', 'R')
    electric_charge = parameter('Electric Charge', 'A s', 'q')

    # Magnetic related parameters
    magnetic_field = parameter('Magnetic Field', 'T', 'B')

    standard_quantities = {
        'length': length,
        'width': width,
        'height': height,
        'diameter': diameter,
        'radius': radius,
        'amplitude': amplitude,
        'thickness': thickness,
        'area': area,
        'volume': volume,
        'mass': mass,
        'velocity': velocity,
        'speed': speed,
        'angular_velocity': angular_velocity,
        'acceleration': acceleration,
        'g': g,
        'force': force,
        'momentum': momentum,
        'period': period,
        'frequency': frequency,
        'energy': energy,
        'work': work,
        'potential_energy': potential_energy,
        'kinetic_energy': kinetic_energy,
        'tension': tension,
        'linear_density': linear_density,
        'stress': stress,
        'power': power,
        'density': density,
        'viscosity': viscosity,
        'pressure': pressure,
        'temperature': temperature,
        'heat': heat,
        'wave_length': wave_length,
        'electric_current': electric_current,
        'electric_voltage': electric_voltage,
        'electric_resistance': electric_resistance,
        'electric_charge': electric_charge,
        'magnetic_field': magnetic_field
    }

    parameters = [standard_quantities[p] for p in parameter_list]
    return parameters


def standard_dimensional_analysis(parameters_string: str) -> dict:
    """Performing dimensional analysis on selected standard physical parameters

    Parameters
    ----------
    parameters_string : str
        name of parameters which will be returned as sympy Quantity object; for
        example: 'density viscosity length'

    Returns
    -------
    dict
        dictionary of dimensionless parameters
    """
    parameters = standard_parameters(parameters_string)
    return dimensional_analysis(*parameters)


def solve_from_standard_dimensional_analysis(parameters_string: str,
                                             target_parameter_string: str
                                             ) -> dict:
    """Performing dimensionall analysis on selected standard physical parameters
    and solving results for specific parameter

    Parameters
    ----------
    parameters_string : str
        name of parameters which will be returned as sympy Quantity object; for
        example: 'density viscosity length'
    target_parameter : sympy Quantity object
        selected parameter for solving result from dimensional analysis

    Returns
    -------
    dict
        dictionary of solutions for each dimensionless parameter
    """
    parameters = standard_parameters(parameters_string)
    target_parameter = standard_parameters(target_parameter_string)[0]
    return solve_from_dimensional_analysis(*parameters,
                                           target_parameter=target_parameter)
