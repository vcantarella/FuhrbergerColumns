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
        k_no3, k_no2, k_no3c, k_no2c, K_no3, K_no2, K_pyr, K_c, δpyr, δc, α, k_dec = p
        # unpack the state variables
        no3_ = u[:,1]
        no2_ = u[:,2]
        s_edc = u[:,3] # solid ed
        s_c = u[:,4] # solid c
        ep = u[:,5] # enzyme for pyr
        ec = u[:,6] # enzyme for c

        δno3 = 200
        δno2 = 400
        # Activation function:
        F_no3 = @. no3_ / (K_no3 + no3_)
        F_no2 = @. no2_ / (K_no2 + no2_)
        Fp = @. s_edc / (s_edc + K_pyr)
        Fc = @. s_c / (s_c + K_c)
        Fpno3 = @. F_no3 * Fp
        Fpno2 = @. F_no2 * Fp
        Fpno3c = @. F_no3 * Fp * Fc
        Fpno2c = @. F_no2 * Fp * Fc
        # Activation function:
        l1 = @. (δno3 + δpyr)*k_no3*Fpno3*(δno2 + δpyr)*k_no2*Fpno2
        l2 = @. (δno3 + δc)*k_no3c*Fpno3c*(δno2 + δc)*k_no2c*Fpno2c
        u1 = @. (1-k_dec*(ec + ep*(δpyr*l1/l2)))/(1-α*Fp*l1/l2)
        u1 = @. ifelse(u1>0, u1, 0)
        u2 = @. 1 - u1
        # Enzyme functions
        @. du[:,5] = α*Fp*u1 - k_dec*ep
        @. du[:,6] = α*Fc*u2 - k_dec*ec
        # reactions
        r_no3 = @. k_no3 * Fpno3
        r_no2 = @. k_no2 * Fpno2
        r_no3c = @. k_no3c * Fpno3c
        r_no2c = @. k_no2c * Fpno2c

        # transport
        c_advec = [c_in;u[:,1:nmob]]
        advec = -v .* diff(c_advec, dims=1) ./ dx
        gradc = diff(u[:,1:nmob], dims=1)./dx
        disp = ([gradc; zeros(1, nmob)]-[zeros(1, nmob); gradc]).* De ./ dx
        # du/dt
        @. du[:,1] = advec[:,1] + disp[:,1] - 7*r_no3 - 2*r_no3c
        @. du[:,2] = advec[:,2] + disp[:,2] + 7*r_no3 + 2*r_no3c - 14*r_no2 - 4*r_no2c
        @. du[:,3] = -r_no3 - 3*r_no2
        @. du[:,4] = -r_no3c - 3*r_no2c
    end
    return rhscyber!
end


function create_simplefredsmodel(v, De, dx, c_in, nmob)
    # Defining the reaction rates model
    # define the ode model

    function rhs!(du, u, p, t)
        # unpack the parameters
        k_no3, k_no2, K_no3, = p
        # unpack the state variables
        no3_ = @view u[:,1]
        no2_ = @view u[:,2]

        # Activation function:

        n_rows = size(u, 1)
        # transport
        # Calculate transport terms directly without temporary arrays
        @inbounds for j in 1:nmob
            # First cell (boundary condition)
            du[1,j] = -v * (u[1,j] - c_in[1,j]) / dx
            
            # Calculate dispersion at first cell - only forward gradient
            grad_fwd = (u[2,j] - u[1,j]) / dx
            du[1,j] += De[1,j] * grad_fwd / dx  # Remove the gradient difference
            
            # Interior cells
            for i in 2:n_rows-1
                # Advection
                du[i,j] = -v * (u[i,j] - u[i-1,j]) / dx
                
                # Dispersion
                grad_fwd = (u[i+1,j] - u[i,j]) / dx
                grad_bwd = (u[i,j] - u[i-1,j]) / dx
                du[i,j] += De[1,j] * (grad_fwd - grad_bwd) / dx
            end
            
            # Last cell
            du[n_rows,j] = -v * (u[n_rows,j] - u[n_rows-1,j]) / dx
            grad_bwd = (u[n_rows,j] - u[n_rows-1,j]) / dx
            du[n_rows,j] += De[1,j] * (0.0 - grad_bwd) / dx  # Zero-gradient at boundary
        end

        @inbounds for k in 1:n_rows
            # Calculate each term once to avoid repeated computation
            r_no3 = @. k_no3 * no3_[k] / (K_no3 + no3_[k])
            r_no2 = @. k_no2 * no2_[k] / (K_no2 + no2_[k])
            
            # Update state variables
            du[k,1] -= r_no3
            du[k,2] += r_no3 - r_no2
        end
    end
    return rhs!
end