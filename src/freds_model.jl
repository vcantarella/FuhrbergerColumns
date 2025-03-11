"""
Model based on the implementation by Marc Franz
"""


"""
    create_fredsmodel(v, De, dx, c_in, nmob)

Creates an ODE right-hand side (RHS) function to solve the coupled reactive-transport problem considering the Fred's model.

# Arguments
- `v::Number`: The advection velocity (m/s)
- `De::Vector`: The dispersion coefficients for each species (m²/s)
- `dx::Number`: The grid spacing (m)
- `c_in::Vector`: The inflow boundary concentrations
- `nmob::Int`: The number of mobile species

# Returns
- A function representing the ODE RHS for the coupled reactive-transport problem

# Example
```julia
rhs = create_fredsmodel(1.0, [0.1, 0.2], 0.01, [1.0, 0.0], 2, 0.4, 2650.0)
```
"""
function create_fredsmodel(v, De, dx, c_in, nmob)
    # Defining the reaction rates model
    # define the ode model

    function rhs!(du, u, p, t)
        # unpack the parameters
        k_no3, k_no2, k_no3c, k_no2c, K_no3, K_no2, K_pyr, K_c, c_t, st = p
        # unpack the state variables
        no3_ = u[:,1]
        no2_ = u[:,2]
        s_edc = u[:,3] # solid ed
        s_c = u[:,4] # solid c

        # Activation function:
        θ = @. 1 / (exp((c_t-s_c)/st/s_c)+1)

        # transport
        c_advec = [c_in;u[:,1:nmob]]
        advec = -v .* diff(c_advec, dims=1) ./ dx
        gradc = diff(u[:,1:nmob], dims=1)./dx
        disp = ([gradc; zeros(1, nmob)]-[zeros(1, nmob); gradc]).* De ./ dx
        # reaction
        r_no3 = @. k_no3 * no3_ / (K_no3 + no3_) * s_edc / (s_edc + K_pyr) * (1-θ)
        r_no2 = @. k_no2 * no2_ / (K_no2 + no2_) * s_edc / (s_edc + K_pyr) * (1-θ)
        r_no3c = @. k_no3c * no3_ / (K_no3 + no3_) * s_c / (s_c + K_c)
        r_no2c = @. k_no2c * no2_ / (K_no2 + no2_) * s_c / (s_c + K_c)
        
        @. du[:,1] = advec[:,1] + disp[:,1] - 7*r_no3 - 2*r_no3c
        @. du[:,2] = advec[:,2] + disp[:,2] + 7*r_no3 + 2*r_no3c - 14*r_no2 - 4*r_no2c
        @. du[:,3] = -r_no3 - 3*r_no2
        @. du[:,4] = -r_no3c - 3*r_no2c
    end
    return rhs!
end



"""
    create_fredsmodel(v, De, dx, c_in, nmob)

Creates an ODE right-hand side (RHS) function to solve the coupled reactive-transport problem considering the Fred's model.

# Arguments
- `v::Number`: The advection velocity (m/s)
- `De::Vector`: The dispersion coefficients for each species (m²/s)
- `dx::Number`: The grid spacing (m)
- `c_in::Vector`: The inflow boundary concentrations
- `nmob::Int`: The number of mobile species

# Returns
- A function representing the ODE RHS for the coupled reactive-transport problem

# Example
```julia
rhs = create_fredsmodel(1.0, [0.1, 0.2], 0.01, [1.0, 0.0], 2, 0.4, 2650.0)
```
"""
function create_cyberneticfredsmodel(v, De, dx, c_in, nmob)
    # Defining the reaction rates model
    # define the ode model

    function rhscyber!(du, u, p, t)
        # unpack the parameters
        k_no3, k_no2, k_no3c, k_no2c, K_no3, K_no2, K_pyr, K_c, δpyr, δc = p
        # unpack the state variables
        no3_ = u[:,1]
        no2_ = u[:,2]
        s_edc = u[:,3] # solid ed
        s_c = u[:,4] # solid c

        # Activation function:


        # transport
        c_advec = [c_in;u[:,1:nmob]]
        advec = -v .* diff(c_advec, dims=1) ./ dx
        gradc = diff(u[:,1:nmob], dims=1)./dx
        disp = ([gradc; zeros(1, nmob)]-[zeros(1, nmob); gradc]).* De ./ dx
        # reaction
        r_no3 = @. k_no3 * no3_ / (K_no3 + no3_) * s_edc / (s_edc + K_pyr)
        r_no2 = @. k_no2 * no2_ / (K_no2 + no2_) * s_edc / (s_edc + K_pyr)
        r_no3c = @. k_no3c * no3_ / (K_no3 + no3_) * s_c / (s_c + K_c)
        r_no2c = @. k_no2c * no2_ / (K_no2 + no2_) * s_c / (s_c + K_c)

        max_r = maximum(hcat(δpyr.*r_no3, δpyr.*r_no2, δc.*r_no3c, δc.*r_no2c), dims=2)
        v_rno3 = @. δpyr*r_no3/max_r
        v_rno2 = @. δpyr*r_no2/max_r
        v_rno3c = @. δc*r_no3c/max_r
        v_rno2c = @. δc*r_no2c/max_r


        @. du[:,1] = advec[:,1] + disp[:,1] - 7*r_no3*v_rno3 - 2*r_no3c*v_rno3c
        @. du[:,2] = advec[:,2] + disp[:,2] + 7*r_no3*v_rno3 + 2*r_no3c*v_rno3c - 14*r_no2*v_rno2 - 4*r_no2c*v_rno2c
        @. du[:,3] = -r_no3*v_rno3 - 3*r_no2*v_rno2
        @. du[:,4] = -r_no3c*v_rno3c - 3*r_no2c*v_rno2c
    end
    return rhscyber!
end