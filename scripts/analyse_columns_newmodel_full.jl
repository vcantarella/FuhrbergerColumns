using DrWatson
using CairoMakie
using DataFrames, XLSX, Statistics
using Dates
using StaticArrays
using OrdinaryDiffEq
using Optimization
using SciMLSensitivity
using SparseConnectivityTracer
using Symbolics
using LinearSolve
# load analytical data
file_path = datadir("exp_raw", "no2_analyticall_results.xlsx")
sheet_name = "sampling"
data = XLSX.readtable(file_path, sheet_name)
df = DataFrame(data)
# start-time and end-time columns are of type Any because some of the values are in DateTiem and some in Time.
# Assume the Time happened on 17.06.2025 and convert them to DateTime
Hour(df[5,"start time"])
start_time = df[!,"start time"]
end_time = df[!,"end time"]
# Start time for the t0 of the experiment for columns 1 and 2
t0_1_2 = DateTime(2025, 05, 12, 20, 0, 0)
c_no3_1_2 = 0.6 # concentration of NO3- in the input solution for columns 1 and 2 [mM]
c_lac_1_2 = 0.0 # concentration of lactate in the input solution for columns 1 and 2 [mM]
# Start time for the t0 of the experiment for columns 3 and 4
t0_3_4 = DateTime(2025, 05, 12, 20, 14, 0)
t0_arr = [ifelse(col == 1 || col == 2, t0_1_2, t0_3_4) for col in df[!,"column"]]
# convert start_time and end_time to seconds since t0_dt
start_time_s = Dates.Second.(start_time - t0_arr)
end_time_s = Dates.Second.(end_time - t0_arr)

# import transport parameters
transp_params = load("data/optimized_tracer_params.jld2")
tracer_params = transp_params["tracer_params"]
# Generate the transport dataset for each column
# Q is the flow rate in μL/min for each column
# And we calculate the v and αₗ*v that are dependent on the current Q
struct TData{QT, T}
    column::Int64
    v::QT # flow velocity in m/s
    De::QT # longitudinal dispersion coefficient in m2/s (without the eff. diffusion)
    end_times::T
end
qvec = TData[]
D = 3.5*1e-2 #cm to m diameter of the column
A = π * D^2 / 4 # Cross-sectional area
for i in 1:4
    # only the values where there are no missing values in the flow rate
    ϕ = tracer_params[i][1]
    αₗ = tracer_params[i][2]
    bool_index = (.!ismissing.(df[:, "Q"])) .& (df[!,"column"] .== i)
    # convert the flow rate to μL/min
    Q = df[bool_index,"Q"] ./ 60 .* 1e-9 # convert from μL/min to m3/s
    v = Q./(ϕ * A) # flow velocity in m/s
    De = v .* αₗ # longitudinal dispersion coefficient in m2/s
    # create a QData object and push it to the vector
    push!(qvec, TData(i, v, De, Dates.value.(end_time_s[bool_index])))
end

struct CinData{QT, T}
    column::Int64
    c_in::QT # inflow concentration in mM
    t_in::T # time of the inflow concentration
end


## Events:
# Start time for the t0 of the experiment for columns 1 and 2
c_no3_1_2 = 0.6 # concentration of NO3- in the input solution for columns 1 and 2 [mM]
c_lac_1_2 = 0.0 # concentration of lactate in the input solution for columns 1 and 2 [mM]
# Start time for the t0 of the experiment for columns 3 and 4
c_no3_3_4 = 0.6 # concentration of NO3- in the input solution for columns 3 and 4 [mM]
c_lac_3_4 = 0.33 # concentration of lactate in the input solution for columns 3 and 4 [mM]
# Moment when we swithed the bag for columns 3 and 4, removing the lactate solution
t_switch_3_4 = DateTime(2025, 05, 16, 16, 35, 0)
c_lac_3_4_switch = 0.0 # concentration of lactate in the input solution for columns 3 and 4 after the switch [mM]
t_s_34_s = Dates.Second(t_switch_3_4 - t0_1_2)
# Moment when we switch the bag for all columns, increasing the NO3- concentration
t_switch_all_1 = DateTime(2025, 05, 23, 19, 00, 0)
c_no3_all_1 = 1.0 # concentration of NO3- in the input solution for all columns after the switch [mM]
t_1_s = Dates.Second(t_switch_all_1 - t0_1_2)
# Moment when the increase the NO3- concentration
t_switch_all_2 = DateTime(2025, 05, 28, 19, 15, 0)
c_no3_all_2 = 2.0 # concentration of NO3- in the input solution for all columns after the switch [mM]
t_2_s = Dates.Second(t_switch_all_2 - t0_1_2)
# Moment when pumps 1 and 2 stopped due to clogging.
t_stop_1 = DateTime(2025,05,28, 19, 50)
t_stop_2 = DateTime(2025,05,28, 20,10)
# moment when the pumps were restarted
t_restart_1_2 = DateTime(2025,05,29, 08, 45)
# Approximate time when all pumps stoppeds due to energie failure in the building
t_stop_all = DateTime(2025,05,29, 12, 30)
t_stop_s = Dates.Second(t_stop_all - t0_1_2)
# Moment when the pumps were restarted after the energy failure
t_restart_all = DateTime(2025,05,29, 18, 16)
t_restart_s = Dates.Second(t_restart_all - t0_1_2)
# Moment when we switch to 4 mM inflow concentration
t_switch_4mM = DateTime(2025, 06, 03, 18, 25, 0)
c_no3_4mM = 4.
t_4_s = Dates.Second(t_switch_4mM - t0_1_2)

c_ins = CinData[]
for i in 1:4
    c_no3 = 0.6e-3
    c_lac = 0.0
    c_so4 = 2.3e-3 # concentration of SO4-2 in the input solution [mM]
    c_fe = 0.0e-3 # concentration of Fe+2 in the input solution [mM]

    if i == 3 || i == 4
        c_lac = 0.33e-3 # concentration of NO3- in the input solution for columns 3 and 4 [mM]
    end
    cins = [[c_no3, 1e-12, c_so4, c_fe, c_lac, c_no3, 0.0],]
    t0s = []
    if i == 3 || i == 4
        # For columns 3 and 4, we have a switch in the lactate concentration
        push!(cins, [c_no3, 1e-12, c_so4, c_fe, c_lac_3_4_switch, c_no3, 0.0])
        push!(t0s, t_s_34_s) # convert days to seconds
        push!(cins, [c_no3_all_1*1e-3, 1e-12, c_so4, c_fe, c_lac_3_4_switch, c_no3_all_1*1e-3, 0.0])
        push!(t0s, t_1_s) # convert days to seconds
        push!(cins, [c_no3_all_2*1e-3, 1e-12, c_so4, c_fe, c_lac_3_4_switch, c_no3_all_2*1e-3, 0.0])
        push!(t0s, t_2_s) # convert days to seconds
        push!(cins, [c_no3_4mM*1e-3, 1e-12, c_so4, c_fe, c_lac_3_4_switch, c_no3_4mM*1e-3, 0.0])
        push!(t0s, t_4_s) # convert days to seconds
    else
        # For columns 1 and 2, we have a switch in the NO3-
        push!(cins, [c_no3_all_1*1e-3, 1e-12, c_so4, c_fe, c_lac, c_no3_all_1*1e-3, 0.0])
        push!(t0s, t_1_s) # convert days to seconds
        push!(cins, [c_no3_all_2*1e-3, 1e-12, c_so4, c_fe, c_lac, c_no3_all_2*1e-3, 0.0])
        push!(t0s, t_2_s) # convert days to seconds
        push!(cins, [c_no3_4mM*1e-3, 1e-12, c_so4, c_fe, c_lac, c_no3_4mM*1e-3, 0.0])
        push!(t0s, t_4_s) # convert days to seconds
    end
    t0s = Dates.value.(t0s) # convert to seconds
    # Add the initial concentration at t0
    cindata = CinData(i, cins, t0s)
    push!(c_ins, cindata)
end

    
# Diffusion coefficients for the different (mobile) components
Deff = [
    1.0e-9, # NO3-
    1.0e-9, # NO2-
    1.0e-9, # SO4-2
    1.0e-9, # Fe+2
    1.0e-9, # Lactate
    1.0e-9, # NO3- (tracer)
]

gs = @SVector [160.16,437.13,-150]
function reactive_transport_builder(flow_data::TData, cin_data::CinData, Deff, dx, ϕ, ρₛ, gs)
    vs = SVector{length(flow_data.v)}(flow_data.v)
    v = mean(vs) # average velocity
    Des = SVector{length(flow_data.De)}(flow_data.De)
    DD = mean(Des) # average dispersion coefficient
    end_times = SVector{length(flow_data.end_times)}(flow_data.end_times)
    c_ins = SVector{length(cin_data.c_in)}(cin_data.c_in)
    t_ins = SVector{length(cin_data.t_in)}(cin_data.t_in)
    De_work = Vector{Float64}(undef, length(Deff))
    De_v2 = Deff .+ DD
    function dynamic_transport(t, De, Deff=Deff)
        for i in eachindex(vs)
            if t <= end_times[i]
                v = vs[i]
                De_i = Des[i]
                De .= Deff .+ De_i
                return v   # This should be inside the if block
            end
        end
        # Add a fallback return outside the loop
        De .= Deff .+ Des[end]  # This allocates a new array every call!
        return vs[end]
    end
    function dynamic_c_in(t)
        for i in eachindex(t_ins)
            if t <= t_ins[i]
                return c_ins[i]
            end
        end
        return c_ins[end]  # Return the last concentration if no match found
    end
    function expo_fun(gs, G₀, GcdG₀, r_no3b, r_no2b, r_so4b)
        # Calculate the exponent for the activation function
        expo = 1 - (gs[1] / G₀ * r_no3b + gs[2] / G₀ * r_no2b + gs[3] / G₀ * r_so4b) - GcdG₀ * (2*r_no3b + 3*r_no2b + 8*r_so4b)
        return expo
    end
    function build_rhs(dynamic_transport, dynamic_c_in, expo_fun, dx, gs,De=De_work, mult_term = (1-ϕ)/ϕ*ρₛ)
        function rhs!(du, u, p, t)
            @inline v = dynamic_transport(t,De)
            # unpack the state variables
            no3_ = @view u[:,1]
            no2_ = @view u[:,2]
            so4_ = @view u[:,3]
            fe2_ = @view u[:,4]
            lac_ = @view u[:,5]
            no3t_ = @view u[:,6]  # NO3- tracer
            cd = @view u[:,7]  # NO3- tracer concentration
            b = @view u[:,8]  # biomass active fraction (immobile))
            #b_in = @view u[:,8]  # biomass inactive fraction (immobile)

            # unpack the parameters
            μ₁, μ₂, αₑ, k_dec, Ya1, Ya2, Yd1, Yd2, Ka1, Ka2, Kd, Kb = p

            n_rows = size(u, 1)
            @inline c_in = dynamic_c_in(t)
            nmob = size(c_ins, 1)  # Number of mobile components
            # transport
            # Calculate transport terms directly without temporary arrays
            @inbounds for j in 1:nmob
                # First cell (boundary condition)
                du[1,j] = -v * (u[1,j] - c_in[j]) / dx
                
                # Calculate dispersion at first cell - only forward gradient
                grad_fwd = (u[2,j] - u[1,j]) / dx
                du[1,j] += De[j] * grad_fwd / dx  # Remove the gradient difference
                
                # Interior cells
                for i in 2:n_rows-1
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
            @inbounds for k in 1:n_rows
                # Calculate each term once to avoid repeated computation
                # γ = αₑ/(b[k]+Kb)*(μ₁/Yd1*no3_[k]/(Ka1 + no3_[k]) + μ₂/Yd2*no2_[k]/(Ka2 + no2_[k]))^(-1)
                # γ = ifelse(γ ≥ 0.9, .9, γ)  # Limit γ to 1
                cd_r = cd[k] / (cd[k] + Kd)
                ca1_r = no3_[k] / (no3_[k] + Ka1)
                ca2_r = no2_[k] / (no2_[k] + Ka2)
                r_no3b = μ₁/Ya1 * ca1_r * cd_r
                r_no2b = μ₂/Ya2 * ca2_r * cd_r
                r_b = ((μ₁ * ca1_r + 
                    μ₂ * ca2_r) * cd_r - k_dec) * b[k]  # Decay term
                r_cd = αₑ * b[k] / (b[k] + Kb) -  μ₁/Yd1*ca1_r * cd_r - μ₂/Yd2*ca2_r * cd_r
                    # Update state variables
                du[k,1] -= r_no3b*b[k]
                du[k,2] += (r_no3b - r_no2b)*b[k]
                # du[k,3] -= mult_term*r_so4b*b[k]
                # du[k,4] += mult_term*r_no2b*b
                # du[k,5] -= mult_term*r_no2b*b

                du[k,7] += r_cd  # biomass active fraction
                du[k,8] = r_b  # biomass inactive fraction
                # du[k,8] += r_ina - r_act  # biomass inactive fraction
            end
        end
    return rhs!
    end
    return build_rhs(dynamic_transport, dynamic_c_in, expo_fun, dx, gs)
end

dx = 0.0005 # Spatial step size
L = 0.08 #m (8 cm)  # Spatial locations
x = range(0+dx/2, stop=L-dx/2, step=dx)  # Spatial locations
rhs! = reactive_transport_builder(qvec[1], c_ins[1], Deff, dx, tracer_params[1][1], 2.65, gs)
γa = 1/25 # stoichiometric coefficient of e-acceptor in the anabolic rctieaction
γc1 = 4/2 # stoichiometric coefficient of e-acceptor in the catabolic rea
γc2 = 4/3 # stoichiometric coefficient of e-acceptor in the catabolic rea
Yd = 0.3 # yield of biomass from the catabolic reaction
p0 = [
    1e-5, # μ₁ (growth rate for NO3-)
    5e-6, # μ₂ (growth rate for NO2-)
    8e-10, # αₑ (effective diffusion coefficient)
    1.15e-6, # k_dec (decay rate)
    1/((1/Yd-1)*γc1),
    1/((1/Yd-1)*γc2),
    Yd,
    Yd,
    5e-4,
    5e-4,
    1e-6,
    5e-5
]

u0 = zeros(length(x), 8) # 5 mobile components + 2 immobile components (active and inactive biomass)
# u0[:,1] .= 1e-12 # Initial concentration of NO3-
# u0[:,2] .= 1e-12 # Initial concentration of NO2-
u0[:,8] .= 1e-7 # Initial concentration of inactive biomass
du0 = copy(zeros(size(u0))) # Initialize the derivative array

rhs!(du0, u0, p0, 0.0) # Calculate the initial derivative
# @time rhs!(du0, u0, p0, 0.0) # Benchmark the RHS function
# using BenchmarkTools
# @benchmark rhs!($du0, $u0, $p0, 0.0) # Benchmark the RHS function
# @code_warntype rhs!(du0, u0, p0, 0.0) # Profile the RHS function

du0
tspan = (0.0, 32*24*60*60) # 32 days in seconds
old_prob = ODEProblem(rhs!, u0, tspan, p0)
# Symbolics for sparsity detection
detector = TracerSparsityDetector()
jac_sparsity2 = ADTypes.jacobian_sparsity((du, u) -> rhs!(du, u, p0, 100000),
    du0, u0, detector) # add the sparsity pattern to speed up the solution
jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> rhs!(du, u, p0, 100000),
    du0, u0) # add the sparsity pattern to speed up the
jac_sparsity == jac_sparsity2
fixed_rhs! = ODEFunction(rhs!, jac_prototype=jac_sparsity2)
fastprob = ODEProblem(fixed_rhs!, u0, tspan, p0)

# defining points to stop
tstops = vcat(qvec[1].end_times, c_ins[1].t_in)
sort!(tstops) # sort the times
sol = solve(fastprob, FBDF(nlsolve=NLNewton(max_iter = 100)), abstol = 1e-8, reltol = 1e-8,
    maxiters = 100000,
    tstops = tstops,
    d_discontinuities = tstops,
    )
sol.t
# Check model outflow:
no3_out = [sol.u[i][end, 1] for i in eachindex(sol.t)]
no2_out = [sol.u[i][end, 2] for i in eachindex(sol.t)]
so4_out = [sol.u[i][end, 3] for i in eachindex(sol.t)]
lac_out = [sol.u[i][end, 5] for i in eachindex(sol.t)]
tracer_out = [sol.u[i][end, 6] for i in eachindex(sol.t)]

# Plot results
fig = Figure()
ax = Axis(fig[1, 1], title = "Outflow concentrations",
    xlabel = "Time (days)", ylabel = "Concentration (M)")
plot_t = sol.t ./ (24*60*60) # convert seconds to days
lines!(ax, plot_t, no3_out, label = "NO3- outflow", color = :blue)
lines!(ax, plot_t, tracer_out, label = "NO3- tracer outflow", color = :blue, linestyle = :dash)
lines!(ax, plot_t, no2_out, label = "NO2- outflow", color = :orange)
lines!(ax, plot_t, so4_out, label = "SO4-2 outflow", color = :green)
lines!(ax, plot_t, lac_out, label = "Lactate outflow", color = :red)
axislegend(ax, position = :lt)
fig

# Plot alive vs inactive biomass
fig2 = Figure()
ax2 = Axis(fig2[1, 1], title = "Biomass concentrations",
    xlabel = "length (m)", ylabel = "Concentration (M)")
plot_x = x ./ 0.01 # convert cm to m
lines!(ax2, plot_x, sol.u[end][:, 7], label = "Active biomass", color = :blue)
#lines!(ax2, plot_x, sol.u[end][:, 7], label = "Inactive biomass", color = :orange)
axislegend(ax2, position = :rt)
fig2
#Dead volumes per column
dvs = Dict(1 => 35,
           2 => 29,
           3 => 34,
           4 => 30) # in cm
dvs_t0 = Dict(1 => 66,
              2 => 40+57,
              3 => 41+38,
              4 => 41+38.5) # in cm
tube_diam = 0.152 # cm
# Dead volume in cm3
dv = Dict(1 => dvs[1] * π * (tube_diam/2)^2,
          2 => dvs[2] * π * (tube_diam/2)^2,
          3 => dvs[3] * π * (tube_diam/2)^2,
          4 => dvs[4] * π * (tube_diam/2)^2)
dv_t0 = Dict(1 => dvs_t0[1] * π * (tube_diam/2)^2,
             2 => dvs_t0[2] * π * (tube_diam/2)^2,
             3 => dvs_t0[3] * π * (tube_diam/2)^2,
             4 => dvs_t0[4] * π * (tube_diam/2)^2)
t0s = Dict(1 => t0_dt + Dates.Second(floor(Int64, dv_t0[1] / Q[1])),
            2 => t0_dt + Dates.Second(floor(Int64, dv_t0[2] / Q[2])),
            3 => t0_dt + Dates.Second(floor(Int64, dv_t0[3] / Q[3])),
            4 => t0_dt + Dates.Second(floor(Int64, dv_t0[4] / Q[4])))
# calculate the start and end
column = convert.(Int64,df[!,"column"])
start_times_dict = Dict(1 => [Dates.value(Dates.Second(t - t0s[1])) for t in start_time[column .== 1]],
                        2 => [Dates.value(Dates.Second(t - t0s[2])) for t in start_time[column .== 2]],
                        3 => [Dates.value(Dates.Second(t - t0s[3])) for t in start_time[column .== 3]],
                        4 => [Dates.value(Dates.Second(t - t0s[4])) for t in start_time[column .== 4]])
end_times_dict = Dict(1 => [Dates.value(Dates.Second(t - t0s[1])) for t in end_time[column .== 1]],
                      2 => [Dates.value(Dates.Second(t - t0s[2])) for t in end_time[column .== 2]],
                      3 => [Dates.value(Dates.Second(t - t0s[3])) for t in end_time[column .== 3]],
                      4 => [Dates.value(Dates.Second(t - t0s[4])) for t in end_time[column .== 4]])
# save the flow rate data to a file
Qs = Dict(1 => Q[column .== 1],
          2 => Q[column .== 2],
          3 => Q[column .== 3],
          4 => Q[column .== 4])
avg_times_dict = Dict(1 => (start_times_dict[1] .+ end_times_dict[1]) ./ 2 .- dv[1] ./ Qs[1],
    2 => (start_times_dict[2] .+ end_times_dict[2]) ./ 2 .- dv[2] ./ Qs[2],
    3 => (start_times_dict[3] .+ end_times_dict[3]) ./ 2 .- dv[3] ./ Qs[3],
    4 => (start_times_dict[4] .+ end_times_dict[4]) ./ 2 .- dv[4] ./ Qs[4]
)