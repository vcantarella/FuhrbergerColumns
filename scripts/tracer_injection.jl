using SpecialFunctions
using DrWatson
using CairoMakie
using DataFrames, XLSX, Statistics
using Dates


# load analytical data
file_path = datadir("exp_raw", "no2_analyticall_results.xlsx")
sheet_name = "br_standard_curve"
data = XLSX.readtable(file_path, sheet_name)
df = DataFrame(data)
# start-time and end-time columns are of type Any because some of the values are in DateTiem and some in Time.
# Assume the Time happened on 17.06.2025 and convert them to DateTime
Hour(df[5,"start-time"])
start_time = DateTime[]
end_time = DateTime[]
for i in eachindex(df[!,"start-time"])
    if typeof(df[i, "start-time"]) == DateTime
        push!(start_time, df[i, "start-time"])
    else #type = Time
        push!(start_time, DateTime(2025, 06, 17, Dates.value(Hour(df[i, "start-time"])),
            Dates.value(Minute(df[i, "start-time"]))))
    end
    if typeof(df[i, "end-time"]) == DateTime
        push!(end_time, df[i, "end-time"])
    else #type = Time
        push!(end_time, DateTime(2025, 06, 17, Dates.value(Hour(df[i, "end-time"])),
            Dates.value(Minute(df[i, "end-time"]))))
    end
end
t0_dt = DateTime(2025, 06, 16, 22, 48) # start time of the experiment
# convert start_time and end_time to seconds since t0_dt
start_time_s = [Dates.Second(t - t0_dt) for t in start_time]
end_time_s = [Dates.Second(t - t0_dt) for t in end_time]
Q = convert.(Float64, df[!,"Q"]) # flow rate in μL/min
Q = Q .* 1e-3 / 60 # convert to cm3/s

# save the flow rate data to a file
save("data/br_flow_rate.jld2", Dict("Q" => Q,
    "start_times" => Dates.value.(start_time_s),
    "end_times" => Dates.value.(end_time_s),
    "column" => convert.(Int64, df[!,"column"])))

struct FlowRate{TQ<:AbstractVector, TTime<:AbstractVector}
    Q::TQ
    start_times::TTime
    end_times::TTime
end

function (f::FlowRate)(t)
    @inbounds for i in eachindex(f.start_times)
        if t >= f.start_times[i] && t <= f.end_times[i]
            return f.Q[i]
        elseif t > f.end_times[i] && t < f.start_times[i+1]
            return (f.Q[i]+ f.Q[i+1])/ 2
        end
    end
    return zero(eltype(f.Q))
end

function calculate_flow_rate(Q, start_time_s, end_time_s)
    # Q is a vector of flow rates in cm3/s
    # start_time_s and end_time_s are vectors of start and end times in seconds since t0_dt
    # return a function that takes a time t and returns the flow rate at that time
    @assert length(Q) == length(start_time_s) == length(end_time_s)
    @assert all(start_time_s .< end_time_s) # start times must be before end times
    @assert all(Q .>= 0) # flow rates must be non-negative
    return FlowRate(Q, start_time_s, end_time_s)
end
# create a FlowRate object per column
column = convert.(Int64,df[!,"column"])
flow_rates = Dict{Int64, FlowRate}()
for col in 1:4
    # get the flow rates for the column
    Q_col = Q[column .== col]
    start_time_s_col = Dates.value.(start_time_s[column .== col])
    end_time_s_col = Dates.value.(end_time_s[column .== col])
    # create a FlowRate object for the column
    flow_rates[col] = calculate_flow_rate(Q_col, start_time_s_col, end_time_s_col)
end
Br = df[!,"Br-"] # concentration in mM
## calculate the dead times to correct the start and end times.
# <To calculate the dead times we use the dead volumes and the flow rates
#Dead volumes per column
dvs = Dict(1 => 66+35,
           2 => 40+57+29,
           3 => 41+38+34,
           4 => 41+38.5+30) # in cm
tube_diam = 0.152 # cm
# Dead volume in cm3
dv = Dict(1 => dvs[1] * π * (tube_diam/2)^2,
          2 => dvs[2] * π * (tube_diam/2)^2,
          3 => dvs[3] * π * (tube_diam/2)^2,
          4 => dvs[4] * π * (tube_diam/2)^2)
# Dead time in seconds
avg_time = @. convert(Float64,Dates.value(start_time_s + (end_time_s - start_time_s) ./ 2))

for i in eachindex(avg_time)
    avg_time[i] -= dv[column[i]] / flow_rates[column[i]](avg_time[i]) # subtract the dead time
end

# Ok now we are ready to create the datasets per column
dfs = Dict{Int64, DataFrame}()
missing_indices = findall(ismissing, Br)
Br = Br[1:end .!= missing_indices]
Br[Br .>= 1] .= 1.0
avg_time = avg_time[1:end .!= missing_indices]
column = column[1:end .!= missing_indices]
for col in 1:4
    # create a DataFrame for the column. Make types concrete and skip missing Br values
    
    avg_time_col = avg_time[column .== col]
    # create a DataFrame with the average time and the Br values for the column
    dfs[col] = DataFrame(
        time = avg_time[column .== col],
        Br = Br[column .== col],
    )
end

# Now we can plot the data
fig = Figure()
ax = Axis(fig[1, 1],
    xlabel = "Time [s]",
    ylabel = "Br⁻ [mM]",
    title = "Bromide breakthrough curves")
for col in 1:4
    # get the data for the column
    df = dfs[col]
    # plot the data
    #lines!(ax, df.time, df.Br, label = "Column $col")
    scatter!(ax, df.time, df.Br, label = "Column $col")
end
# add a legend
axislegend(ax, position = :lt, framevisible = false)
fig
# save the figure
save("plots/br_breakthrough_data.png", fig)

dfs.keys # check the keys of the dictionary
#change the keys to strings
dfs = Dict{String, DataFrame}(string(k) => v for (k, v) in pairs(dfs))
# save the datasets to a file
save("data/br_breakthrough_data.jld2", dfs)

"""
constant_injection(cr, x, t, c0, c_in, v, Dl)

Calculates the concentration profile of a solute undergoing
advection-dispersion transport in a porous 1D domain with constant
 in a 1D domain with constant injection
Original reference: (Ogata & Banks, 1961): 
Appelo, C.A.J.; Postma, Dieke. Geochemistry, Groundwater and Pollution (p. 104).

# Arguments
- `cr::Matrix`: A 2D array to store the concentration profile of a solute. 
    first dimension is time and second dimension is space.
- `x::Vector`: A 1D array of spatial locations. (Where to sample the concentration)
- `t::Vector`: A 1D array of time locations. (When to sample the concentration)
- `c0::Real`: The concentratio at x=0 (inflow concentration).
- `c_in::Real`: The initial concentration in the column (t=0).
- `v::Real`: The velocity of the solute.
- `Dl::Real`: The longitudinal dispersion coefficient.
# Returns
    nothing, the results are stored in the `cr` array. 
"""
function constant_injection!(
     cr::Matrix,
     x::Vector,
     t::Vector,
     c0::Real,
     c_in::Real,
     v::Real,
     Dl::Real,
     )
    
    for i in eachindex(x)
        cr[:, i] .= c_in .+ (c0 - c_in) / 2 .* erfc.((x[i] .- v .* t)
         ./ (2 .* sqrt.(Dl .* t)))
    end
    return nothing
end

Q = 21 # μL/min - imprecise
avg_q = Q*1e-3 / 60 # cm3/s
d = 3.5 # cm
A = π * (d/2)^2 # cm2
L = 8 # cm
x = [L/100,]
t = collect(0:360:(36*3600))
cr = zeros(length(t), length(x))

ϕ₀ = 0.37 # porosity
L/100/(avg_q/A/ϕ₀/100)
αₗ = 1e-3 # m
D = 2e-9 * ϕ₀ # m2/s
q = avg_q/A/100 # m/s
v = q/ϕ₀ # m/s
dv0 = 100*π*0.152^2/4 #cm3
p = [ϕ₀, αₗ]
# function cost(p)
#     q_red, αₗ = p
#     v = q*q_red/ϕ₀
#     # dt = dv/avg_q # in s 
#     De = αₗ*v .+ D
#     c_in = 440
#     c0 = maximum(ec)
#     cr = zeros(length(t_model), length(x))
#     # t_model_fun = ifelse.(t_model.-dt .< 0, repeat([0.], length(t_model)), t_model.-dt)
#     constant_injection!(cr, x, t_model, c0, c_in, v, De)
#     return sum(abs2,(ec-cr[:, end]))
# end

# cost(p)

# using PRIMA
# xl = [0.4, 1e-8]
# up = [1, 1e-1]
# res = bobyqa(cost, p, xl=xl, xu=up, rhobeg = 1e-1)
# p = res[1]
ϕ, αₗ = p
v = q/ϕ
De = αₗ*v .+ D
c_in = 1e-9
c0 = 1e-3
constant_injection!(cr, x, t, c0, c_in, v, De)
fig = Figure()
ax = Axis(fig[1, 1],
    xlabel = "Time [hr]",
    ylabel = "Br⁻ [mol/L]",
    title = "Bromide breakthrough fit")
lines!(ax, t./3600, cr[:, end], color = :blue)
# scatter!(ax, t, ec, color = :red)
fig
# save("plots/porosity_fit.png", fig)