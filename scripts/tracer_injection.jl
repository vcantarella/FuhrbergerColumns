using SpecialFunctions
using CairoMakie
using DataFrames, XLSX, Statistics
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