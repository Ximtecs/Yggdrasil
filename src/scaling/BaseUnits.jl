"""
    BaseUnits

A mutable struct representing the base physical units in either the SI or CGS system.
"""
mutable struct BaseUnits
    system::String  # The unit system, either 'SI', 'CGS', or 'HL'
    # Physical constants
    eps_0::Float64  # Permittivity of free space (units depend on the system)
    mu_0::Float64   # Permeability of free space (units depend on the system)
    c::Float64      # Speed of light (m/s for SI, cm/s for CGS)
    h_p::Float64    # Planck constant (J·s for SI, erg·s for CGS)
    k_B::Float64    # Boltzmann constant (J/K for SI, erg/K for CGS)
    m_e::Float64    # Electron mass (kg for SI, g for CGS)
    m_p::Float64    # Proton mass (kg for SI, g for CGS)
    m_u::Float64    # Atomic unit mass (kg for SI, g for CGS)
    e::Float64      # Elementary charge (C for SI, StatC for CGS)
    G::Float64      # Gravitational constant (m^3·kg^-1·s^-2 for SI, cm^3·g^-1·s^-2 for CGS)
    eV::Float64     # Electron volt (J for SI, erg for CGS)
    #-------
    eps_0_scaling::Float64 
    mu_0_scaling::Float64 
    charge_scaling::Float64 
    m_e_scaling::Float64 
    """
        BaseUnits(system::String)
    Constructor for BaseUnits. Initializes physical constants based on the specified system ('SI' or 'CGS').
    This method checks if the provided system is valid ('SI' or 'CGS', or 'hL) and initializes the struct with the appropriate physical constants.
    """
    function BaseUnits(system::String; ϵ_s::Float64=1.0, μ_s::Float64=1.0, me_s ::Float64=1.0, q_s::Float64=1.0)
        new_system = system in ["SI", "CGS", "HL"] ? system : throw(Exception("Unit system must be either 'SI', 'CGS', or 'HL"))
        obj = new(new_system)
        set_base_physical_values(obj)  # Set the physical constants based on the system
        #----------- set default scaling to 1.0 (no scaling) ----------------
        scale_base_units(obj, ϵ_s, μ_s, me_s, q_s)
        #---------------------------------------------------------------------
        return obj
    end
end

"""
    set_base_physical_values(obj::BaseUnits)

Sets the base physical values based on the system specified in the object.
"""
function set_base_physical_values(obj::BaseUnits)
    if obj.system == "SI"
        set_base_physical_values_SI(obj)
    elseif obj.system == "CGS"
        set_base_physical_values_CGS(obj)
    elseif obj.system == "HL"
        set_base_physical_values_HL(obj)
    else
        throw(Exception("Invalid unit system"))
    end
end


"""
    set_base_physical_values_SI(obj::BaseUnits)

Sets the base physical values for the SI (International System of Units) unit system.
This function initializes the physical constants in the BaseUnits object with their respective values in the SI system.
"""
function set_base_physical_values_SI(obj::BaseUnits)
    # Permittivity of free space in F/m (farads per meter)
    obj.eps_0 = 8.8541878128e-12
    # Permeability of free space in H/m (henrys per meter)
    obj.mu_0  = 1.25663706212e-6
    # Speed of light in m/s (meters per second), calculated from eps_0 and mu_0
    obj.c = 1 / sqrt(obj.eps_0 * obj.mu_0)
    # Planck constant in J·s (joule seconds)
    obj.h_p = 6.62607015e-34
    # Boltzmann constant in J/K (joules per kelvin)
    obj.k_B = 1.38064852e-23
    # Electron mass in kg (kilograms)
    obj.m_e = 9.10938356e-31
    # Proton mass in kg (kilograms)
    obj.m_p = 1.672621898e-27
    # Atomic unit mass in kg (kilograms)
    obj.m_u = 1.66053906660e-27
    # Elementary charge in C (coulombs)
    obj.e   = 1.6021766208e-19
    # Gravitational constant in m^3·kg^-1·s^-2 (cubic meters per kilogram per square second)
    obj.G   = 6.67430e-11
    # Electron Volt in J (joules)
    obj.eV  = 1.6021772e-19
end

"""
    set_base_physical_values_CGS(obj::BaseUnits)

Sets the base physical values for the CGS (Centimeter-Gram-Second) unit system.
This function initializes the physical constants in the BaseUnits object with their respective values in the CGS system.
"""
function set_base_physical_values_CGS(obj::BaseUnits)
    # Permittivity of free space in s^2·cm^-2 (seconds squared per centimeter squared)
    obj.eps_0 = 1.1126500560536184e-21 /  (4 * pi ) 
    # Permeability of free space, unitless in CGS
    obj.mu_0  = 4 * pi #1.0
    # Speed of light in cm/s (centimeters per second), calculated from eps_0 and mu_0
    obj.c = 1 / sqrt(obj.eps_0 * obj.mu_0)
    # Planck constant in erg·s (erg seconds)
    obj.h_p = 6.62607015e-27
    # Boltzmann constant in erg/K (ergs per kelvin)
    obj.k_B = 1.38064852e-16
    # Electron mass in g (grams)
    obj.m_e = 9.10938356e-28
    # Proton mass in g (grams)
    obj.m_p = 1.672621898e-24
    # Atomic unit mass in g (grams)
    obj.m_u = 1.66053906660e-24
    # Elementary charge in StatC (statcoulombs, CGS unit of charge)
    obj.e   = 4.80320427e-10
    # Gravitational constant in cm^3·g^-1·s^-2 (cubic centimeters per gram per square second)
    obj.G   = 6.67430e-8
    # Electron Volt in erg (ergs)
    obj.eV  = 1.6021772e-12
end

"""
    set_base_physical_values_CGS(obj::BaseUnits)

Sets the base physical values for the CGS (Centimeter-Gram-Second) unit system.
This function initializes the physical constants in the BaseUnits object with their respective values in the CGS system.
"""
function set_base_physical_values_HL(obj::BaseUnits)
    # Permittivity of free space in s^2·cm^-2 (seconds squared per centimeter squared)
    obj.eps_0 = 1.1126500560536184e-21 
    # Permeability of free space, unitless in CGS
    obj.mu_0  = 1.0
    # Speed of light in cm/s (centimeters per second), calculated from eps_0 and mu_0
    obj.c = 1 / sqrt(obj.eps_0 * obj.mu_0)
    # Planck constant in erg·s (erg seconds)
    obj.h_p = 6.62607015e-27
    # Boltzmann constant in erg/K (ergs per kelvin)
    obj.k_B = 1.38064852e-16
    # Electron mass in g (grams)
    obj.m_e = 9.10938356e-28
    # Proton mass in g (grams)
    obj.m_p = 1.672621898e-24
    # Atomic unit mass in g (grams)
    obj.m_u = 1.66053906660e-24
    # Elementary charge in StatC (statcoulombs, CGS unit of charge)
    obj.e   = 4.80320427e-10  * sqrt(4 * pi)
    # Gravitational constant in cm^3·g^-1·s^-2 (cubic centimeters per gram per square second)
    obj.G   = 6.67430e-8
    # Electron Volt in erg (ergs)
    obj.eV  = 1.6021772e-12
end


"""
    print_all(obj::BaseUnits)

Prints all the base physical values of the `BaseUnits` object.
This function checks the system (SI or CGS) of the `BaseUnits` object and calls the appropriate printing function to display the values.
"""
function print_all(obj::BaseUnits)
    # Check the unit system of the object and call the corresponding print function
    if obj.system == "SI"
        print_all_SI(obj)  # Call the function to print SI units
    elseif obj.system == "CGS"
        print_all_CGS(obj)  # Call the function to print CGS units
    end
end

"""
    print_all_SI(obj::BaseUnits)

Prints all the SI base physical values from the `BaseUnits` object.
"""
function print_all_SI(obj::BaseUnits)
    @printf("SI base physical values:\n")
    @printf("     eps_0 = %.4e [F/m] (Permittivity of free space)\n", obj.eps_0)
    @printf("      mu_0 = %.4e [H/m] (Permeability of free space)\n", obj.mu_0)
    @printf("         c = %.4e [m/s] (Speed of light)\n", obj.c)
    @printf("       h_p = %.4e [J·s] (Planck constant)\n", obj.h_p)
    @printf("       k_B = %.4e [J/K] (Boltzmann constant)\n", obj.k_B)
    @printf("       m_e = %.4e [kg] (Electron mass)\n", obj.m_e)
    @printf("       m_p = %.4e [kg] (Proton mass)\n", obj.m_p)
    @printf("       m_u = %.4e [kg] (Atomic unit mass)\n", obj.m_u)
    @printf("         e = %.4e [C] (Elementary charge)\n", obj.e)
    @printf("         G = %.4e [m^3·kg^-1·s^-2] (Gravitational constant)\n", obj.G)
    @printf("        eV = %.4e [J] (Electron Volt)\n", obj.eV)
end

"""
    print_all_CGS(obj::BaseUnits)

Prints all the CGS base physical values from the `BaseUnits` object.
"""
function print_all_CGS(obj::BaseUnits)
    @printf("CGS base physical values:\n")
    @printf("     eps_0 = %.4e [s^2·cm^-2] (Permittivity of free space)\n", obj.eps_0)
    @printf("      mu_0 = %.4e [unitless] (Permeability of free space)\n", obj.mu_0)
    @printf("         c = %.4e [cm/s] (Speed of light)\n", obj.c)
    @printf("       h_p = %.4e [erg·s] (Planck constant)\n", obj.h_p)
    @printf("       k_B = %.4e [erg/K] (Boltzmann constant)\n", obj.k_B)
    @printf("       m_e = %.4e [g] (Electron mass)\n", obj.m_e)
    @printf("       m_p = %.4e [g] (Proton mass)\n", obj.m_p)
    @printf("       m_u = %.4e [g] (Atomic unit mass)\n", obj.m_u)
    @printf("         e = %.4e [StatC] (Elementary charge)\n", obj.e)
    @printf("         G = %.4e [cm^3·g^-1·s^-2] (Gravitational constant)\n", obj.G)
    @printf("        eV = %.4e [erg] (Electron Volt)\n", obj.eV)
end


"""
    print_all_HL(obj::BaseUnits)

Prints all the Heaviside-Lorentz base physical values from the `BaseUnits` object.
"""
function print_all_HL(obj::BaseUnits)
    @printf("CGS base physical values:\n")
    @printf("     eps_0 = %.4e [s^2·cm^-2] (Permittivity of free space)\n", obj.eps_0)
    @printf("      mu_0 = %.4e [unitless] (Permeability of free space)\n", obj.mu_0)
    @printf("         c = %.4e [cm/s] (Speed of light)\n", obj.c)
    @printf("       h_p = %.4e [erg·s] (Planck constant)\n", obj.h_p)
    @printf("       k_B = %.4e [erg/K] (Boltzmann constant)\n", obj.k_B)
    @printf("       m_e = %.4e [g] (Electron mass)\n", obj.m_e)
    @printf("       m_p = %.4e [g] (Proton mass)\n", obj.m_p)
    @printf("       m_u = %.4e [g] (Atomic unit mass)\n", obj.m_u)
    @printf("         e = %.4e [StatC] (Elementary charge)\n", obj.e)
    @printf("         G = %.4e [cm^3·g^-1·s^-2] (Gravitational constant)\n", obj.G)
    @printf("        eV = %.4e [erg] (Electron Volt)\n", obj.eV)
end

"""
    scale_base_units(units::BaseUnits, eps_0_scaling::Float64, mu_0_scaling::Float64, m_e_scaling::Float64, charge_scaling::Float64)

Scales specific base units by the given scaling factors.
"""
function scale_base_units(units::BaseUnits, eps_0_scaling::Float64, mu_0_scaling::Float64, m_e_scaling::Float64, charge_scaling::Float64)
    units.eps_0 *= eps_0_scaling
    units.mu_0 *= mu_0_scaling
    units.e *= charge_scaling
    units.m_e *= m_e_scaling
    units.eps_0_scaling = eps_0_scaling
    units.mu_0_scaling = mu_0_scaling
    units.charge_scaling = charge_scaling
    units.m_e_scaling = m_e_scaling
    units.c = 1 / sqrt(units.eps_0 * units.mu_0)
end
