#----------- strcut for code units needed to reconstuct real units ----
struct code_units
    c::Float64
    B::Float64
    p::Float64
    rho::Float64
    debye::Float64
    skin_depth::Float64
end
#------------------------------------------------

#-------- struct with the single real value needed to reconstruct all real units ----
struct known_real_units
    number_density ::Float64
end
#------------------------------------------------


#-------- struct with the single real value needed to reconstruct all real units ----
struct real_unit_input
    number_density ::Float64
    mass_density :: Float64
    B_field ::Float64
    T :: Float64
end
#------------------------------------------------


#-------- strcut for fudging physical valeus -----------
struct ratios_and_scaling
    electron_proton_ratio :: Float64
    eps_0_scaling         :: Float64
    q_charge_scaling      :: Float64
end
#------------------------------------------------


function find_real_units_CGS(code_units, known_real_units, ratios_and_scaling)
    
    base = BaseUnits("CGS")
    #------------- Scale base units ----------------
        #------ calculate electron mass scaing based on ratio of electron to proton mass -----
    electron_mass_scaling = (base.m_p / base.m_e) / ratios_and_scaling.electron_proton_ratio
        #-------------------------------------------------------------------------------------
    charge_scaling = ratios_and_scaling.q_charge_scaling
    eps_0_scaling = ratios_and_scaling.eps_0_scaling
    mu_0_scaling = 1e0    
    scale_base_units(base, eps_0_scaling, mu_0_scaling, electron_mass_scaling, charge_scaling)
    #-----------------------------------------------


    #------------ physical constants depending on scaling ----------------
    mu = ( (base.m_e + base.m_p)  / 2. ) /  base.m_u# assuming only protons and electrons
    e = base.e # may be scaled
    c_real = base.c
    #------------------------------------------------------------------------
    
    k_B = base.k_B
    
    #--------- calculate real mass density ---------------
    number_density = known_real_units.number_density
    rho_real = number_density * 2 * mu * base.m_u
    #------------------------------------------------
    
    
    
    #------------ Load code units ----------------
    p_code          = code_units.p
    rho_code        = code_units.rho
    B_code          = code_units.B
    c_code          = code_units.c
    debye_code      = code_units.debye
    skin_depth_code = code_units.skin_depth
    #------------------------------------------------
    

    #------------- Calculate real values and scaling factors ----------------
    T = (p_code / rho_code ) * ( (mu * base.m_u) / k_B) * (c_real / c_code)^2
    p_real = ( rho_real / (mu * base.m_u) )  * k_B * T 
    mass_density_scale = rho_real / rho_code
    B_real = B_code * sqrt(p_real / p_code)



    if (debye_code > 0.0 && skin_depth_code > 0.0)
        error("Both debye_code and skin_depth cannot be larger than 1.")
    end

    if (debye_code > 0.0 )
        length_scale = 1. / (debye_code * e) * sqrt(k_B * T / (4 * pi * number_density )) # This is for Gaussian units where mu_0 = 4 * pi
    elseif (skin_depth_code > 0.0)
        skin_depth_real = c_real / (e * sqrt(4 * pi * number_density / base.m_p))
        length_scale = skin_depth_real / skin_depth_code
    else
        error("Either debye_code or skin_depth must be larger than 0.")
    end 



    time_scale = length_scale * c_code / c_real 
    #----------------------------------------------------------------
    

    #------------- assume equal electron and proton temperature, and no temperature scale ----
    T_e = T
    T_p = T
    temperature_scale = 1e0 # don't scale temperature for now - can be added later
    #-----------------------------------------------------------------------------------------

    #---------- Assume commen length is the same as length scaling ---  
    length = length_scale
    #----------------------------------------------------------------


    #---------- create a new scaling object containing all the scaling factors ----
    scaling = ScalingCGS(base,
                number_density, length, T, T_e, T_p, B_real,
                length_scale, mass_density_scale, time_scale, temperature_scale)
    #--------------------------------------------------------------------------------
    return scaling
end

function find_real_units_HL(code_units, known_real_units, ratios_and_scaling)
    
    base = BaseUnits("HL")
    #------------- Scale base units ----------------
        #------ calculate electron mass scaing based on ratio of electron to proton mass -----
    electron_mass_scaling = (base.m_p / base.m_e) / ratios_and_scaling.electron_proton_ratio
        #-------------------------------------------------------------------------------------
    charge_scaling = ratios_and_scaling.q_charge_scaling
    eps_0_scaling = ratios_and_scaling.eps_0_scaling
    mu_0_scaling = 1e0    
    scale_base_units(base, eps_0_scaling, mu_0_scaling, electron_mass_scaling, charge_scaling)
    #-----------------------------------------------


    #------------ physical constants depending on scaling ----------------
    mu = ( (base.m_e + base.m_p)  / 2. ) /  base.m_u# assuming only protons and electrons
    e = base.e # may be scaled
    c_real = base.c
    #------------------------------------------------------------------------
    
    k_B = base.k_B
    
    #--------- calculate real mass density ---------------
    number_density = known_real_units.number_density
    rho_real = number_density * 2 * mu * base.m_u
    #------------------------------------------------
    
    
    
    #------------ Load code units ----------------
    p_code     = code_units.p
    rho_code   = code_units.rho
    B_code     = code_units.B
    c_code     = code_units.c
    debye_code = code_units.debye
    skin_depth_code = code_units.skin_depth
    #------------------------------------------------
    

    #------------- Calculate real values and scaling factors ----------------
    T = (p_code / rho_code ) * ( (mu * base.m_u) / k_B) * (c_real / c_code)^2
    p_real = ( rho_real / (mu * base.m_u) )  * k_B * T 
    mass_density_scale = rho_real / rho_code
    B_real = B_code * sqrt(p_real / p_code)


    #obj.e_real * sqrt(obj.number_density_real / obj.m_p_real)
    skin_depth_real = c_real / (e * sqrt(number_density / base.m_p))




    if (debye_code > 0.0 && skin_depth > 0.0)
        error("Both debye_code and skin_depth cannot be larger than 1.")
    end

    if (debye_code > 0.0 )
        length_scale = 1. / (debye_code * e) * sqrt(k_B * T / (number_density )) # This is for Gaussian units where mu_0 = 4 * pi
    elseif (skin_depth_code > 0.0)
        skin_depth_real = c_real / (e * sqrt(number_density / base.m_p))
        length_scale = skin_depth_real / skin_depth_code
    else
        error("Either debye_code or skin_depth must be larger than 1.")
    end 

    time_scale = length_scale * c_code / c_real 
    #----------------------------------------------------------------
    

    #------------- assume equal electron and proton temperature, and no temperature scale ----
    T_e = T
    T_p = T
    temperature_scale = 1e0 # don't scale temperature for now - can be added later
    #-----------------------------------------------------------------------------------------

    #---------- Assume commen length is the same as length scaling ---  
    length = length_scale
    #----------------------------------------------------------------


    #---------- create a new scaling object containing all the scaling factors ----
    scaling = ScalingHL(base,
                number_density, length, T, T_e, T_p, B_real,
                length_scale, mass_density_scale, time_scale, temperature_scale)
    #--------------------------------------------------------------------------------
    return scaling
end
