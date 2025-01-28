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
        k_no3, k_no2, K_no3, K_no2, c_t, st, k_act = p
        # unpack the state variables
        no3_ = u[:,1]
        no2_ = u[:,2]
        s_edc = u[:,3] # solid ed
        b = u[:,4] # biomass

        # Activation function:
        θ = @. 1 / (exp((c_t-no3_)/st/no3_)+1)
        bmax = 1

        # transport
        c_advec = [c_in;u[:,1:nmob]]
        advec = -v .* diff(c_advec, dims=1) ./ dx
        gradc = diff(u[:,1:nmob], dims=1)./dx
        disp = ([gradc; zeros(1, nmob)]-[zeros(1, nmob); gradc]).* De ./ dx
        # reaction
        r_no3 = @. k_no3 * no3_ * (K_no3 + no3_) * b
        r_no2 = @. k_no2 * no2_ * (K_no2 + no2_) * b
        
        du[:,1] .= advec[:,1] .+ disp[:,1] .- r_no3
        du[:,2] .= advec[:,2] .+ disp[:,2] .+ r_no3 .- r_no2
        du[:,3] .= .-r_no3.*2 .- r_no2.*3
        @. du[:,4] = θ*k_act*(1 -b/bmax)- (1-θ)*k_act
    end
    return rhs!
end