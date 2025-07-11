

function create_simplefredsmodel(v, De, dx, c_in, nmob)
    # Defining the reaction rates model where we dynamically update the boundary conditions
    # of flow rate and inflow concentration.

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