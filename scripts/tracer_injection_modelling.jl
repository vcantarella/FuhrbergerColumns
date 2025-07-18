using OrdinaryDiffEq
using Optimization
using SciMLSensitivity
using SparseConnectivityTracer
using LinearSolve
using JLD2
using CairoMakie
using PRIMA

function create_tracer_rhs_q_callback(D, c_in, dx)
    function tracer_primitive!(du, u, p, t)
        n_rows = size(u, 1)
        ϕ = p[1]
        αₗ = p[2]
        q = p[3]  # Flow rate as a function of time
        # q_l = Float64(q.(t))
        v = q / ϕ  # Velocity based on flow rate and porosity
        De = D + αₗ * v  # Dispersion coefficient
        
        # Calculate transport terms directly without temporary arrays
        # First cell (boundary condition)
        du[1] = -v * (u[1] - c_in) / dx
        
        # Calculate dispersion at first cell
        grad_fwd = (u[2] - u[1]) / dx
        du[1] += De * (grad_fwd - zero(eltype(u))) / dx
        
        # Interior cells
        for i in 2:n_rows-1
            # Advection
            du[i] = -v * (u[i] - u[i-1]) / dx
            
            # Dispersion
            grad_fwd = (u[i+1] - u[i]) / dx
            grad_bwd = (u[i] - u[i-1]) / dx
            du[i] += De * (grad_fwd - grad_bwd) / dx
        end
        
        # Last cell
        du[n_rows] = -v * (u[n_rows] - u[n_rows-1]) / dx
        grad_bwd = (u[n_rows] - u[n_rows-1]) / dx
        du[n_rows] += De * (zero(eltype(u)) - grad_bwd) / dx  # Zero-gradient at boundary
    end
    return tracer_primitive!
end


function create_tracer_rhs_q_function(D, c_in, dx, q)
    function tracer_primitive!(du, u, p, t)
        n_rows = size(u, 1)
        ϕ = p[1]
        αₗ = p[2]
        #q = p[3]  # Flow rate as a function of time
        # q_l = Float64(q.(t))
        v = q(t) / ϕ  # Velocity based on flow rate and porosity
        De = D + αₗ * v  # Dispersion coefficient
        
        # Calculate transport terms directly without temporary arrays
        # First cell (boundary condition)
        du[1] = -v * (u[1] - c_in) / dx
        
        # Calculate dispersion at first cell
        grad_fwd = (u[2] - u[1]) / dx
        du[1] += De * (grad_fwd - zero(eltype(u))) / dx
        
        # Interior cells
        for i in 2:n_rows-1
            # Advection
            du[i] = -v * (u[i] - u[i-1]) / dx
            
            # Dispersion
            grad_fwd = (u[i+1] - u[i]) / dx
            grad_bwd = (u[i] - u[i-1]) / dx
            du[i] += De * (grad_fwd - grad_bwd) / dx
        end
        
        # Last cell
        du[n_rows] = -v * (u[n_rows] - u[n_rows-1]) / dx
        grad_bwd = (u[n_rows] - u[n_rows-1]) / dx
        du[n_rows] += De * (zero(eltype(u)) - grad_bwd) / dx  # Zero-gradient at boundary
    end
    return tracer_primitive!
end
# Define parameters
dx = 0.0001 # Spatial step size
L = 0.08 #m (8 cm)  # Spatial locations
D = 3.5*1e-2 #cm to m
A = π * D^2 / 4 # Cross-sectional area
x = range(0+dx/2, stop=L-dx/2, step=dx)  # Spatial locations
t = range(0.1, stop=72000, length=100)  # Time locations
c_in = 1e-3
ϕ = 0.3
αₗ = 8e-5
De = 1e-9
rhs! = create_tracer_rhs_q_callback(De, c_in, dx)

# load the data and plot the results
dic = load("data/br_breakthrough_data.jld2")
Br_dic = dic["Br"]
time_exp = dic["avg_time"]

ps = Dict(
    1 => [ϕ, αₗ],
    2 => [ϕ, αₗ],
    3 => [ϕ, αₗ],
    4 => [ϕ, αₗ],
)

large_fig = Figure()
large_ax = Axis(large_fig[1, 1],
    xlabel = "Time [hr]",
    ylabel = "Br⁻ [mol/L]",
    title = "Bromide breakthrough fit for all columns")
colors = [:crimson, :steelblue, :forestgreen, :darkorange]

flow_rate_ds = load("data/br_flow_rate.jld2")
for column in [1, 2, 3, 4]
    global ϕ, αₗ  # Declare as global to access/modify global variables
    #column = 1
    Q_1 = flow_rate_ds["Qs"][column]
    start_times = convert.(Float64, flow_rate_ds["start_times_dict"][column])
    end_times = convert.(Float64, flow_rate_ds["end_times_dict"][column])
    Q_1 = Q_1/1e6
    q = Q_1 / A  # Initial flow rate in m^3/s
    # v_est = q / ϕ  # Initial velocity in m/s
    p_tracer = [ϕ, αₗ, q[1]]
    u0 = zeros(length(x))
    du0 = copy(u0)


    tspan = (0.0, maximum(end_times))
    function switch_callback_condition!(out, u, t, integrator, end_times)
        for i in eachindex(end_times)
            if i == 1
                out[i] = t <= end_times[i] ? 0 : 1 # t < end_times[1]
            else
                out[i] = end_times[i-1] <= t <= end_times[i] ? 0 : 1 # t < end_times[2]
            end
        end
    end

    # function q_f(t)
    #     for i in eachindex(end_times)
    #         t <= end_times[i] && return q[i]
    #     end
    #     return 0.0 # Return zero if t is beyond the last end time
    # end
    switch_callback_condition!(out, u, t, integrator) = switch_callback_condition!(out, u, t, integrator, end_times)
    out = zeros(length(end_times))
    switch_callback_condition!(out, u0, 100_000, nothing) # Initialize the condition
    affect!(integrator, event_idx, q_a) = integrator.p[3] = q_a[event_idx]
    affect!(integrator, event_idx) = affect!(integrator, event_idx, q)
    len = length(q)
    cb = VectorContinuousCallback(switch_callback_condition!, affect!, len)
    # rhs! = create_tracer_rhs_q_function(D, c_in, dx, q_f)
    old_prob = ODEProblem(rhs!, u0, tspan, p_tracer)
    # Symbolics for sparsity detection
    detector = TracerSparsityDetector()
    jac_sparsity2 = ADTypes.jacobian_sparsity((du, u) -> rhs!(du, u, p_tracer, 1),
        du0, u0, detector) # add the sparsity pattern to speed up the solution
    fixed_rhs! = ODEFunction(rhs!, jac_prototype=jac_sparsity2)
    fastprob = ODEProblem(fixed_rhs!, u0, tspan, p_tracer)
    sol = solve(fastprob, QNDF(linsolve=KrylovJL_GMRES()), callback = cb, abstol = 1e-8, reltol = 1e-8,
        maxiters = 10000,
        tstops = end_times,
        )
    br = Br_dic[column]*1e-3
    time = time_exp[column]
    # fig = Figure()
    # ax = Axis(fig[1, 1],
    #     xlabel = "Time [hr]",
    #     ylabel = "Br⁻ [mol/L]",
    #     title = "Bromide breakthrough fit")
    # lines!(ax, sol.t ./ 3600, [sol.u[i][end] for i in 1:length(sol.u)], color = :red, label = "Model Fit")
    # scatter!(ax, time ./ 3600, br, color = :blue, label = "Experimental Data")
    # fig

    ## Optimizing the parameters

    function cost(p)
        model_p = [p;[q[1]]]   
        prob = remake(fastprob, p = model_p)
        sol = solve(prob, QNDF(linsolve= KrylovJL_GMRES()), callback = cb, abstol = 1e-8, reltol = 1e-8,
            saveat = time,
            maxiters = 10000,
            tstops = end_times,
            )
        # find the index that match time_exp with sol.t
        valid_idx = findall(x-> x in time, sol.t)
        result = [sol.u[i][end] for i in valid_idx]
        return sum(abs2, (br - result))
    end

    cost([ϕ, αₗ])


    xl = [0.1, 1e-4]
    up = [0.5, 1e-2]
    res = bobyqa(cost, [ϕ, αₗ], xl=xl, xu=up)
    p = res[1]
    ps[column] = p
    ϕ, αₗ = p
    p_new = [ϕ, αₗ, q[1]]
    # Update the problem with the new parameters
    fastprob = remake(fastprob, p = p_new)
    sol = solve(fastprob, QNDF(linsolve= KrylovJL_GMRES()), callback = cb, abstol = 1e-10, reltol = 1e-10,
        maxiters = 10000,
        tstops = end_times,
        )
    # fig = Figure()
    # ax = Axis(fig[1, 1],
    #     xlabel = "Time [hr]",
    #     ylabel = "Br⁻ [mol/L]",
    #     title = "Bromide breakthrough fit")
    lines!(large_ax, sol.t ./ 3600, [sol.u[i][end] for i in 1:length(sol.u)], color = colors[column], label = "Column $column Model Fit")
    scatter!(large_ax, time ./ 3600, br, color = colors[column], label = "Column $column Data")
end
axislegend(large_ax, position = :lt, "Columns")
large_fig
save("plots/br_breakthrough_fit.png", large_fig)

# Save the parameters
save("data/optimized_tracer_params.jld2", "tracer_params", ps)