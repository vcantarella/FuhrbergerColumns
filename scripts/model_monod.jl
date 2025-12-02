
function monod_model(v, De, dx, nmob)
    # Defining the reaction rates model where we dynamically update the boundary conditions
    # of flow rate and inflow concentration.
    function rhs!(du, u, p, t)
        # unpack the parameters
        μ_m, k_dec, K_no3, Bmax, Ya = p
        # unpack the state variables
        no3_ = @view u[:,1]
        b = @view u[:,2]

        # unpack the state variables
        c_in = @view u[1, :] # first row is the inflow concentration

        # Activation function:
        n_rows = size(u, 1)
        # transport
        # Calculate transport terms directly without temporary arrays
        @inbounds for j in 1:nmob
                # First cell (boundary condition)
                du[2,j] = -v * (u[2,j] - c_in[j]) / dx
                
                # Calculate dispersion at first cell - only forward gradient
                grad_fwd = (u[3,j] - u[2,j]) / dx
                du[2,j] += De[j] * grad_fwd / dx  # Remove the gradient difference

                # Interior cells
                for i in 3:n_rows-1
                    # Advection
                    du[i,j] = -v * (u[i,j] - u[i-1,j]) / dx
                    
                    # Dispersion
                    grad_fwd = (u[i+1,j] - u[i,j]) / dx
                    grad_bwd = (u[i,j] - u[i-1,j]) / dx
                    du[i,j] += De[j] * (grad_fwd - grad_bwd) / dx
                end
                
                # Last cell
                du[n_rows,j] = -v * (u[n_rows,j] - u[n_rows-1,j]) / dx
                grad_bwd = (u[n_rows,j] - u[n_rows-1,j]) / dx
                du[n_rows,j] += De[j] * (0.0 - grad_bwd) / dx  # Zero-gradient at boundary
            end

        @inbounds for k in 2:n_rows
            μ = μ_m * no3_[k] / (no3_[k] + K_no3) * b[k]
            # Update state variables
            du[k,1] -= μ/Ya
            du[k,2] = μ*(1-b[k]/Bmax) - k_dec*b[k]
        end
        #make sure du[1, :] = 0
        du[1, :] .= 0.0
    end
    return rhs!
end
function zero_order_model(v, De, dx, nmob)
    function rhs!(du, u, p, t)
            # unpack the parameters
            r_s, = p
            n_rows = size(u, 1) # number of spatial rows excluding the inflow row
            c_in = @view u[1, :] # first row is the inflow concentration
            # transport
            # Calculate transport terms directly without temporary arrays
            @inbounds for j in 1:nmob
                # First cell (boundary condition)
                du[2,j] = -v * (u[2,j] - c_in[j]) / dx
                
                # Calculate dispersion at first cell - only forward gradient
                grad_fwd = (u[3,j] - u[2,j]) / dx
                du[2,j] += De[j] * grad_fwd / dx  # Remove the gradient difference

                # Interior cells
                for i in 3:n_rows-1
                    # Advection
                    du[i,j] = -v * (u[i,j] - u[i-1,j]) / dx
                    
                    # Dispersion
                    grad_fwd = (u[i+1,j] - u[i,j]) / dx
                    grad_bwd = (u[i,j] - u[i-1,j]) / dx
                    du[i,j] += De[j] * (grad_fwd - grad_bwd) / dx
                end
                
                # Last cell
                du[n_rows,j] = -v * (u[n_rows,j] - u[n_rows-1,j]) / dx
                grad_bwd = (u[n_rows,j] - u[n_rows-1,j]) / dx
                du[n_rows,j] += De[j] * (0.0 - grad_bwd) / dx  # Zero-gradient at boundary
            end
            @inbounds for k in 2:n_rows
                r_no3 = ifelse(u[k, 1] > 0, r_s, 0.0) # constant term for NO3-
                       # Update state variables
                du[k,1] -= r_no3
            end
            #make sure du[1, :] = 0
            du[1, :] .= 0.0
        end
    return rhs!
end