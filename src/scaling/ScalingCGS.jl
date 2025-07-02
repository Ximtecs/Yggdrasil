
"""
    ScalingCGS

A struct for handling scaling in the CGS (Centimeter-Gram-Second) unit system.
"""
mutable struct ScalingCGS
    system::String  # Unit system ('CGS' in this case)

    # Input physical parameters
    number_density_real::Float64       # Number density [cm^-3]
    length_real::Float64               # Characteristic length [cm]
    temperature_real::Float64          # Temperature [K]
    temperature_e_real::Float64        # Electron Temperature [K]
    temperature_p_real::Float64        # Proton Temperature [K]
    B_flux_real::Float64               # Magnetic flux density [Gauss]

    # Scaling factors
    length_scaling::Float64            # Length scaling factor
    mass_density_scaling::Float64      # Mass density scaling factor
    time_scaling::Float64              # Time scaling factor
    temperature_scaling :: Float64     # Temperature scaling factor

    # Base physical constants (from BaseUnits)
    base_units::BaseUnits

    # Derived real values
    mu_real::Float64                   # Mean molecular weight [unitless]
    mass_density_real::Float64         # Mass density [g/cm^3]
    time_real::Float64                 # Time [s]
    u_real::Float64                    # Characteristic velocity [cm/s]
    v_a_real::Float64                  # Alfven velocity [cm/s]
    v_thermal_e_real::Float64          # Electron thermal velocity [cm/s]
    v_thermal_p_real::Float64          # Proton thermal velocity [cm/s]
    electron_gyro_freq_real::Float64   # Electron gyro frequency [rad/s]
    proton_gyro_freq_real::Float64     # Proton gyro frequency [rad/s]
    electron_plasma_freq_real::Float64 # Electron plasma frequency [rad/s]
    proton_plasma_freq_real::Float64   # Proton plasma frequency [rad/s]
    skin_depth_real::Float64           # Electron skin depth [cm]
    skin_depth_p_real::Float64         # Proton skin depth [cm]
    debye_len_real::Float64            # Debye length [cm]
    inter_dist_real::Float64           # Interparticle distance [cm]
    e_gyro_radiues_real::Float64       # Electron gyro radius [cm]
    p_gyro_radiues_real::Float64       # Proton gyro radius [cm]
    mass_real::Float64                 # Characteristic mass [g]
    energy_real::Float64               # Characteristic energy [erg]
    energy_per_mass_real::Float64      # Energy per mass [erg/g]
    energy_per_volume_real::Float64    # Energy per volume [erg/cm^3]
    pressure_real::Float64             # Pressure [Ba]


    eps_0_real::Float64  # Permittivity of free space (units depend on the system)
    mu_0_real::Float64   # Permeability of free space (units depend on the system)
    c_real::Float64      # Speed of light (m/s for SI, cm/s for CGS)
    h_p_real::Float64    # Planck constant (J·s for SI, erg·s for CGS)
    k_B_real::Float64    # Boltzmann constant (J/K for SI, erg/K for CGS)
    m_e_real::Float64    # Electron mass (kg for SI, g for CGS)
    m_p_real::Float64    # Proton mass (kg for SI, g for CGS)
    m_u_real::Float64    # Atomic unit mass (kg for SI, g for CGS)
    e_real::Float64      # Elementary charge (C for SI, StatC for CGS)
    G_real::Float64      # Gravitational constant (m^3·kg^-1·s^-2 for SI, cm^3·g^-1·s^-2 for CGS)
    eV_real::Float64     # Electron volt (J for SI, erg for CGS)


    # Additional fields for scaling factors and other properties can be added here

    #scaling values 
    mu_scaling :: Float64
    number_density_scaling::Float64
    u_scaling::Float64
    v_a_scaling::Float64
    v_thermal_e_scaling::Float64
    v_thermal_p_scaling::Float64
    electron_gyro_freq_scaling::Float64
    proton_gyro_freq_scaling::Float64
    electron_plasma_freq_scaling::Float64
    proton_plasma_freq_scaling::Float64
    skin_depth_scaling::Float64
    skin_depth_p_scaling::Float64
    debye_len_scaling::Float64
    inter_dist_scaling::Float64
    e_gyro_radiues_scaling::Float64
    p_gyro_radiues_scaling::Float64
    mass_scaling::Float64
    B_flux_scaling :: Float64
    #temperature_scaling::Float64
    temperature_e_scaling::Float64
    temperature_p_scaling::Float64
    energy_scaling::Float64
    energy_per_mass_scaling::Float64
    energy_per_volume_scaling::Float64
    pressure_scaling::Float64


    eps_0_scaling::Float64  
    mu_0_scaling::Float64  
    c_scaling::Float64       
    h_p_scaling::Float64    
    k_B_scaling::Float64     
    m_e_scaling::Float64    
    m_p_scaling::Float64    
    m_u_scaling::Float64    
    e_scaling::Float64      
    G_scaling::Float64      
    eV_scaling::Float64     

    #code values 
    mu_code :: Float64
    mass_density_code::Float64
    number_density_code::Float64
    u_code::Float64
    v_a_code::Float64
    v_thermal_e_code::Float64
    v_thermal_p_code::Float64
    time_code :: Float64
    electron_gyro_freq_code::Float64
    proton_gyro_freq_code::Float64
    electron_plasma_freq_code::Float64
    proton_plasma_freq_code::Float64
    length_code :: Float64
    skin_depth_code::Float64
    skin_depth_p_code::Float64
    debye_len_code::Float64
    inter_dist_code::Float64
    e_gyro_radiues_code::Float64
    p_gyro_radiues_code::Float64
    mass_code::Float64
    B_flux_code :: Float64
    temperature_code::Float64
    temperature_e_code::Float64
    temperature_p_code::Float64
    energy_code::Float64
    energy_per_mass_code::Float64
    energy_per_volume_code::Float64
    pressure_code::Float64


    eps_0_code::Float64  
    mu_0_code::Float64  
    c_code::Float64       
    h_p_code::Float64    
    k_B_code::Float64     
    m_e_code::Float64    
    m_p_code::Float64    
    m_u_code::Float64    
    e_code::Float64      
    G_code::Float64      
    eV_code::Float64     


    #For Macro particle weights
    ds :: Float64
    per_cell :: Float64
    weight_scaled :: Float64
    weight_scaled_mass_proton :: Float64 
    weight_scaled_mass_electron :: Float64 
    weight_scaled_charge :: Float64 
    rho_n_ratio :: Float64 
    
    #Lorentz-Maxwell Unit conversions
    k_E  :: Float64
    k_B  :: Float64 
    k_F  :: Float64 
    k_D  :: Float64
    k_M  :: Float64 
    k_H  :: Float64



    # function ScalingCGS(base_units, number_density=1e15, length=1e8, temperature=1e6, temperature_e=1e6, temperature_p=1e6, B_flux=1e1, length_scale=3e8, mass_density_scale=1e-13, time_scale=1.0, temperature_scale=1e6)
    #     system = base_units.system
    #     obj = new(system, number_density, length, temperature, temperature_e, temperature_p, B_flux, length_scale, mass_density_scale, time_scale,temperature_scale, base_units)
    #     # Additional initialization and calculations for derived values
    #     copy_base_units(obj)
    #     set_derived_real_values(obj)
    #     set_scaling_factors(obj)
    #     calculate_code_values(obj)
    #     set_Maxwell_Lorentz(obj)        
    #     return obj
    # end

    function ScalingCGS(base_units :: BaseUnits, n = 1e15, l=1e8, T=1e6, Te = 1e6, Tp = 1e6, B = 10, l_s = 1e8, ρ_s = 1e-13, t_s  = 100, Temp_s = 1, ds = 0.1, PPC = 128 )
        system = base_units.system
        obj = new(system, n, l, T, Te, Tp, B, l_s, ρ_s, t_s,Temp_s, base_units)
        copy_base_units(obj)
        set_derived_real_values(obj)
        set_scaling_factors(obj)
        calculate_code_values(obj)
        set_Maxwell_Lorentz(obj)     
        set_macro_particle_weights(obj, ds, PPC)
        return obj
    end 


    function ScalingCGS(base_units :: BaseUnits, n = 1e15, l=1e8, T=1e6, Te = 1e6, Tp = 1e6, B = 10, l_s = 1e8, ρ_s = 1e-13, t_s  = 100, Temp_s = 1)
        system = base_units.system
        obj = new(system, n, l, T, Te, Tp, B, l_s, ρ_s, t_s,Temp_s, base_units)
        copy_base_units(obj)
        set_derived_real_values(obj)
        set_scaling_factors(obj)
        calculate_code_values(obj)
        set_Maxwell_Lorentz(obj)     
        return obj
    end 

end

"""
    copy_base_units

Copies the base unit from the BaseUnits struct to the ScalingCGS struct
"""
function copy_base_units(obj::ScalingCGS)
    # Copy base unit values from BaseUnits to ScalingCGS
    obj.eps_0_real = obj.base_units.eps_0   # Permittivity of free space [unitless in CGS]
    obj.mu_0_real  = obj.base_units.mu_0    # Permeability of free space [s^2 cm^-2 in CGS]
    obj.c_real     = obj.base_units.c       # Speed of light [cm s^-1]
    obj.h_p_real   = obj.base_units.h_p     # Planck constant [erg s]
    obj.k_B_real   = obj.base_units.k_B     # Boltzmann constant [erg/K]
    obj.m_e_real   = obj.base_units.m_e     # Electron mass [g]
    obj.m_p_real   = obj.base_units.m_p     # Proton mass [g]
    obj.m_u_real   = obj.base_units.m_u     # Atomic unit mass [g]
    obj.e_real     = obj.base_units.e       # Elementary charge [StatC]
    obj.G_real     = obj.base_units.G       # Gravitational Constant [cm^3 g^-1 s^-2]
    obj.eV_real    = obj.base_units.eV      # Electron Volt [erg]
end


"""
    set_derived_real_values

Calcualtes all real values used in the scaling struct
    Depends on the input parameters: 
        - number_density
        - length
        - temperature
        - temperature_e
        - temperature_p
        - B_flux
"""
function set_derived_real_values(obj::ScalingCGS)
    #Calculate mu_real
    # Note -we calculate this to account for scaled electron masses which may change mu
    obj.mu_real = ( (obj.m_e_real + obj.m_p_real)  / 2. ) /  obj.m_u_real # assuming only protons and electrons

    # Calculate mass density [g/cm^3]
    obj.mass_density_real = obj.number_density_real * 2 * obj.mu_real * obj.m_u_real

    # Calculate time scale [s] - Length scale divided by Alfven velocity
    obj.time_real = obj.length_real / (obj.B_flux_real / sqrt(4 * π * obj.mass_density_real))

    # Calculate characteristic velocity [cm/s]
    obj.u_real = obj.length_real / obj.time_real


    # Calculate Alfven velocity [cm/s] and rescale... 
    obj.v_a_real = obj.B_flux_real / sqrt(4 * π * obj.mass_density_real)
    obj.v_a_real *= obj.c_real / sqrt(obj.v_a_real^2 + obj.c_real^2) # rescaling only changes the actual value if v_a is close to the speed of light
    #TODO - double check scaling does not change anything

    # Calculate electron and proton thermal velocities [cm/s]
    obj.v_thermal_e_real = sqrt(obj.k_B_real * obj.temperature_e_real / obj.m_e_real)
    obj.v_thermal_p_real = sqrt(obj.k_B_real * obj.temperature_p_real / obj.m_p_real)

    # Calculate frequencies [rad/s]
    obj.electron_gyro_freq_real = obj.e_real * obj.B_flux_real / (obj.m_e_real * obj.c_real)
    obj.proton_gyro_freq_real = obj.e_real * obj.B_flux_real / (obj.m_p_real * obj.c_real)
    obj.electron_plasma_freq_real = obj.e_real * sqrt(4 * π * obj.number_density_real / obj.m_e_real)
    obj.proton_plasma_freq_real = obj.e_real * sqrt(4 * π * obj.number_density_real / obj.m_p_real)

    # Calculate lengths [cm]
    obj.skin_depth_real = obj.c_real / obj.electron_plasma_freq_real
    obj.skin_depth_p_real = obj.c_real / obj.proton_plasma_freq_real
    obj.debye_len_real = obj.v_thermal_e_real / obj.electron_plasma_freq_real
    obj.inter_dist_real = 1 / obj.number_density_real^(1/3)
    obj.e_gyro_radiues_real = obj.v_thermal_e_real / obj.electron_gyro_freq_real
    obj.p_gyro_radiues_real = obj.v_thermal_p_real / obj.proton_gyro_freq_real

    # Calculate other derived values
    obj.mass_real = obj.mass_density_real * obj.length_real^3                         # [ g ]
    obj.energy_real = obj.mass_real * obj.u_real^2                                    # [ erg ]
    obj.energy_per_mass_real = obj.energy_real / obj.mass_real                        # [ erg  g^-1]
    obj.energy_per_volume_real = obj.energy_real / obj.length_real^3                  # [ erg  cm^-3]
    #multiply with 2 to account for both protons and electrons
    obj.pressure_real = obj.number_density_real * obj.k_B_real * obj.temperature_real * 2 # [ Ba ]
end

function set_scaling_factors(obj::ScalingCGS)
    # Length scalings
    obj.skin_depth_scaling = obj.length_scaling
    obj.skin_depth_p_scaling = obj.length_scaling
    obj.debye_len_scaling = obj.length_scaling
    obj.inter_dist_scaling = obj.length_scaling
    obj.e_gyro_radiues_scaling = obj.length_scaling
    obj.p_gyro_radiues_scaling = obj.length_scaling

    # Frequency scalings
    obj.electron_gyro_freq_scaling = 1.0 / obj.time_scaling
    obj.proton_gyro_freq_scaling = 1.0 / obj.time_scaling
    obj.electron_plasma_freq_scaling = 1.0 / obj.time_scaling
    obj.proton_plasma_freq_scaling = 1.0 / obj.time_scaling

    # Velocity scalings
    obj.u_scaling = obj.length_scaling / obj.time_scaling
    obj.c_scaling = obj.u_scaling
    obj.v_a_scaling = obj.u_scaling
    obj.v_thermal_e_scaling = obj.u_scaling
    obj.v_thermal_p_scaling = obj.u_scaling

    # Other scalings
    obj.mass_scaling = obj.mass_density_scaling * obj.length_scaling^3
    obj.energy_scaling = obj.mass_scaling * obj.length_scaling^2 / obj.time_scaling^2
    obj.energy_per_mass_scaling = obj.length_scaling^2 / obj.time_scaling^2
    obj.energy_per_volume_scaling = obj.mass_scaling / obj.length_scaling / obj.time_scaling^2
    obj.pressure_scaling = obj.energy_per_volume_scaling

    obj.B_flux_scaling = 1.0 / sqrt(obj.length_scaling) * sqrt(obj.mass_scaling) / obj.time_scaling
    obj.number_density_scaling = 1.0 / obj.length_scaling^3

    obj.temperature_scaling = obj.temperature_scaling
    obj.temperature_e_scaling = obj.temperature_scaling
    obj.temperature_p_scaling = obj.temperature_scaling
    obj.mu_scaling = 1.0

    # Fundamental constants scaling
    obj.G_scaling = obj.length_scaling^3 / obj.mass_scaling / obj.time_scaling^2
    obj.k_B_scaling = obj.energy_scaling / obj.temperature_scaling
    obj.h_p_scaling = obj.energy_scaling * obj.time_scaling
    obj.m_u_scaling = obj.mass_scaling
    obj.m_p_scaling = obj.mass_scaling
    obj.m_e_scaling = obj.mass_scaling
    obj.e_scaling = obj.length_scaling^(3.0/2.0) * sqrt(obj.mass_scaling) / obj.time_scaling
    obj.eps_0_scaling = obj.time_scaling^2 / obj.length_scaling^2
    obj.mu_0_scaling = 1.0  # No scaling
    obj.eV_scaling = obj.energy_scaling
end

function calculate_code_values(obj::ScalingCGS)
    # Fundamental constants code values
    obj.G_code = obj.G_real / obj.G_scaling
    obj.k_B_code = obj.k_B_real / obj.k_B_scaling
    obj.h_p_code = obj.h_p_real / obj.h_p_scaling
    obj.m_u_code = obj.m_u_real / obj.m_u_scaling
    obj.m_p_code = obj.m_p_real / obj.m_p_scaling
    obj.m_e_code = obj.m_e_real / obj.m_e_scaling

    obj.e_code = obj.e_real / obj.e_scaling
    obj.eps_0_code = obj.eps_0_real / obj.eps_0_scaling
    obj.mu_0_code = obj.mu_0_real / obj.mu_0_scaling
    obj.eV_code = obj.eV_real / obj.eV_scaling

    # Length code values
    obj.skin_depth_code = obj.skin_depth_real / obj.skin_depth_scaling
    obj.skin_depth_p_code = obj.skin_depth_p_real / obj.skin_depth_p_scaling
    obj.debye_len_code = obj.debye_len_real / obj.debye_len_scaling
    obj.inter_dist_code = obj.inter_dist_real / obj.inter_dist_scaling
    obj.e_gyro_radiues_code = obj.e_gyro_radiues_real / obj.e_gyro_radiues_scaling
    obj.p_gyro_radiues_code = obj.p_gyro_radiues_real / obj.p_gyro_radiues_scaling

    # Frequency code values
    obj.electron_gyro_freq_code = obj.electron_gyro_freq_real / obj.electron_gyro_freq_scaling
    obj.proton_gyro_freq_code = obj.proton_gyro_freq_real / obj.proton_gyro_freq_scaling
    obj.electron_plasma_freq_code = obj.electron_plasma_freq_real / obj.electron_plasma_freq_scaling
    obj.proton_plasma_freq_code = obj.proton_plasma_freq_real / obj.proton_plasma_freq_scaling

    # Velocity code values
    obj.u_code = obj.u_real / obj.u_scaling
    obj.c_code = obj.c_real / obj.c_scaling
    obj.v_a_code = obj.v_a_real / obj.v_a_scaling
    obj.v_thermal_e_code = obj.v_thermal_e_real / obj.v_thermal_e_scaling
    obj.v_thermal_p_code = obj.v_thermal_p_real / obj.v_thermal_p_scaling

    # Other code values
    obj.mass_code = obj.mass_real / obj.mass_scaling
    obj.energy_code = obj.energy_real / obj.energy_scaling
    obj.energy_per_mass_code = obj.energy_per_mass_real / obj.energy_per_mass_scaling
    obj.energy_per_volume_code = obj.energy_per_volume_real / obj.energy_per_volume_scaling
    obj.pressure_code = obj.pressure_real / obj.pressure_scaling
    obj.B_flux_code = obj.B_flux_real / obj.B_flux_scaling
    obj.number_density_code = obj.number_density_real / obj.number_density_scaling
    obj.temperature_code = obj.temperature_real / obj.temperature_scaling
    obj.temperature_e_code = obj.temperature_e_real / obj.temperature_e_scaling
    obj.temperature_p_code = obj.temperature_p_real / obj.temperature_p_scaling

    obj.mu_code = obj.mu_real / obj.mu_scaling

    obj.mass_density_code = obj.mass_density_real / obj.mass_density_scaling
    obj.length_code = obj.length_real / obj.length_scaling
    obj.time_code = obj.time_real / obj.time_scaling
end



function set_macro_particle_weights(obj::ScalingCGS, ds, per_cell)
    obj.ds = ds
    obj.per_cell = per_cell

    obj.weight_scaled = obj.number_density_code * ds^3 / per_cell

    obj.weight_scaled_mass_proton = obj.weight_scaled * obj.m_p_code
    obj.weight_scaled_mass_electron = obj.weight_scaled * obj.m_e_code
    obj.weight_scaled_charge = obj.weight_scaled * obj.e_code

    obj.rho_n_ratio = 1. / (2. * obj.mu_code * obj.m_u_code)
end


function set_Maxwell_Lorentz(obj:: ScalingCGS)
    obj.k_E = 1. 
    obj.k_B = 1. / obj.c_code
    obj.k_F = 1. / obj.c_code

    obj.k_D = 1. 
    obj.k_M = obj.c_code
    obj.k_H = 1.
end 


function print_Maxwell_Lorentz(obj :: ScalingCGS)
    @printf("\n Maxwell Lorents Factors used in code:\n")
    @printf("%-50s = % .4e \n", " k_E = ", obj.k_E)
    @printf("%-50s = % .4e \n", " k_B = ", obj.k_B)
    @printf("%-50s = % .4e \n", " k_F = ", obj.k_F)
    @printf("%-50s = % .4e \n", " k_D = ", obj.k_D)
    @printf("%-50s = % .4e \n", " k_M = ", obj.k_M)
    @printf("%-50s = % .4e \n", " k_H = ", obj.k_H)
end

function print_Macro_particle_weights(obj :: ScalingCGS)
    @printf("\n Macro particle weights:\n")
    @printf("%-50s = % .4e \n", "                          ds = ", obj.ds           )
    @printf("%-50s = % .4e \n", "                    per_cell = ", obj.per_cell     )
    @printf("%-50s = % .4e \n", "               weight_scaled = ", obj.weight_scaled)
    @printf("%-50s = % .4e \n", " weight_scaled_mass_electron = ", obj.weight_scaled_mass_electron)
    @printf("%-50s = % .4e \n", "   weight_scaled_mass_proton = ", obj.weight_scaled_mass_proton)
    @printf("%-50s = % .4e \n", "        weight_scaled_charge = ", obj.weight_scaled_charge)
    @printf("%-50s = % .4e \n", "                 rho_n_ratio = ", obj.rho_n_ratio)
end

function print_fundamentals(obj::ScalingCGS)
    @printf("\n Fundamental physical constants:\n")
    @printf("%-50s = % .4e % .4e % .4e\n", "        Gravitational Constant [ cm^3 g^-1 s^-2 ]", obj.G_real, obj.G_scaling, obj.G_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "          Boltzmann's Constant [ erg K^-1 ]", obj.k_B_real, obj.k_B_scaling, obj.k_B_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "             Planck's Constant [ erg s ]", obj.h_p_real, obj.h_p_scaling, obj.h_p_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "              Atomic Mass Unit [ g ]", obj.m_u_real, obj.m_u_scaling, obj.m_u_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "                   Proton Mass [ g ]", obj.m_p_real, obj.m_p_scaling, obj.m_p_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "                 Electron Mass [ g ]", obj.m_e_real, obj.m_e_scaling, obj.m_e_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "             Elementary Charge [ StatC ]", obj.e_real, obj.e_scaling, obj.e_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "           Vacuum Permittivity [ s^2 cm^-2 ]", obj.eps_0_real, obj.eps_0_scaling, obj.eps_0_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "           Vacuum Permeability [  ]", obj.mu_0_real, obj.mu_0_scaling, obj.mu_0_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "                Speed of Light [ cm/s ]", obj.c_real, obj.c_scaling, obj.c_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "                 Electron Volt [ erg ]", obj.eV_real, obj.eV_scaling, obj.eV_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "         Mean molecular weight [ ]", obj.mu_real, obj.mu_scaling, obj.mu_code)
end


function print_defining_constants(obj::ScalingCGS)
    @printf("\n Defining Constants:\n")
    @printf("%-50s = % .4e % .4e % .4e\n", "                 Number Desity [ cm^3 ]", obj.number_density_real, obj.number_density_scaling, obj.number_density_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "         Characteristic Length [ cm ]", obj.length_real, obj.length_scaling, obj.length_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "                   Temperature [ K ]", obj.temperature_real, obj.temperature_scaling, obj.temperature_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "          Electron Temperature [ K ]", obj.temperature_e_real, obj.temperature_e_scaling, obj.temperature_e_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "            Proton Temperature [ K ]", obj.temperature_p_real, obj.temperature_p_scaling, obj.temperature_p_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "        Magnetic Flux Strength [ G ]", obj.B_flux_real, obj.B_flux_scaling, obj.B_flux_code)
end


function print_velocities(obj::ScalingCGS)
    @printf("\n Velocities:\n")
    @printf("%-50s = % .4e % .4e % .4e\n", "       Characteristic velocity [ cm s^-1 ]", obj.u_real, obj.u_scaling, obj.u_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "               Alfven velocity [ cm s^-1 ]", obj.v_a_real, obj.v_a_scaling, obj.v_a_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "     electron thermal velocity [ cm s^-1 ]", obj.v_thermal_e_real, obj.v_thermal_e_scaling, obj.v_thermal_e_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "       proton thermal velocity [ cm s^-1 ]", obj.v_thermal_p_real, obj.v_thermal_p_scaling, obj.v_thermal_p_code)
end


function print_frequencies(obj::ScalingCGS)
    @printf("\n Frequencies:\n")
    @printf("%-50s = % .4e % .4e % .4e\n", "     electron gyro frequency   [ rad s^-1 ]", obj.electron_gyro_freq_real, obj.electron_gyro_freq_scaling, obj.electron_gyro_freq_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "       proton gyro frequency   [ rad s^-1 ]", obj.proton_gyro_freq_real, obj.proton_gyro_freq_scaling, obj.proton_gyro_freq_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "     electron plasma frequency [ rad s^-1 ]", obj.electron_plasma_freq_real, obj.electron_plasma_freq_scaling, obj.electron_plasma_freq_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "       proton plasma frequency [ rad s^-1 ]", obj.proton_plasma_freq_real, obj.proton_plasma_freq_scaling, obj.proton_plasma_freq_code)
end



function print_lengths(obj::ScalingCGS)
    @printf("\n Lengths:\n")
    @printf("%-50s = % .4e % .4e % .4e\n", "           electron skin depth [ cm ]", obj.skin_depth_real, obj.skin_depth_scaling, obj.skin_depth_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "             proton skin depth [ cm ]", obj.skin_depth_p_real, obj.skin_depth_p_scaling, obj.skin_depth_p_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "                 debye Length  [ cm ]", obj.debye_len_real, obj.debye_len_scaling, obj.debye_len_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "       inter particle distance [ cm ]", obj.inter_dist_real, obj.inter_dist_scaling, obj.inter_dist_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "          electron gyro radius [ cm ]", obj.e_gyro_radiues_real, obj.e_gyro_radiues_scaling, obj.e_gyro_radiues_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "            proton gyro radius [ cm ]", obj.p_gyro_radiues_real, obj.p_gyro_radiues_scaling, obj.p_gyro_radiues_code)
end

function print_other(obj::ScalingCGS)
    @printf("\n Other:\n")
    @printf("%-50s = % .4e % .4e % .4e\n", "                  Mass density [ g cm^-3 ]", obj.mass_density_real, obj.mass_density_scaling, obj.mass_density_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "           Characteristic mass [ g ]", obj.mass_real, obj.mass_scaling, obj.mass_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "          Characteristic time  [ s ]", obj.time_real, obj.time_scaling, obj.time_code)
    @printf("%-50s = % .4e % .4e % .4e\n", " Characteristic kinetic energy [ erg ]", obj.energy_real, obj.energy_scaling, obj.energy_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "               energy per mass [ erg g^-1]", obj.energy_per_mass_real, obj.energy_per_mass_scaling, obj.energy_per_mass_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "             energy per volume [ erg cm^-3]", obj.energy_per_volume_real, obj.energy_per_volume_scaling, obj.energy_per_volume_code)
    @printf("%-50s = % .4e % .4e % .4e\n", "      Characteristic  pressure [ Ba ]", obj.pressure_real, obj.pressure_scaling, obj.pressure_code)
end

function print_all_CGS(obj::ScalingCGS)
    print_Maxwell_Lorentz(obj)
    print_Macro_particle_weights(obj)
    @printf("\n\t\t\t\t\t\t\tReal      Scaling      Code")
    print_fundamentals(obj)
    print_defining_constants(obj)
    print_velocities(obj)
    print_frequencies(obj)
    print_lengths(obj)
    print_other(obj)
end

function print_basic_info(obj::ScalingCGS)
    println("lenght                 = ", obj.length_real)
    println("number_density         = ", obj.number_density_real)
    println("B                      = ", obj.B_flux_real)
    println("T                      = ", obj.temperature_real)
    println("")
    println("length_scale           = ", obj.length_scaling)
    println("mass_density_scale     = ", obj.mass_density_scaling)
    println("time_scale             = ", obj.time_scaling)
    println("temperature_scale      = ", obj.temperature_scaling)
    println("")
    println("e_mass_scale           = ", obj.base_units.m_e_scaling)
    println("charge_scaling         = ", obj.base_units.charge_scaling)
    println("eps_0_scaling          = ", obj.base_units.eps_0_scaling)
    println("mu_0_scaling           = ", obj.base_units.mu_0_scaling)
end
