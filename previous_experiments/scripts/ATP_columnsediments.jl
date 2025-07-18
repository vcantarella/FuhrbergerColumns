using Statistics
using CairoMakie
using Interpolations
using Colors

# Hydrogeology Group - Dinithi Amarawardana 
# Date 2024.06.18 === 2024: 06:28

# STANDARDS 
# were made using the same artifical ground water used for the
# experiment and it was autoclaved. 100uM stock solution of ATP (Thermo
# Scientific ATP solution) was used to make the standards of atp
# assay was carried out using BacTiter-Glow microbial cell viability
# assay (Promega)

# Sampling - 

# 2.5 g of each was taken in 15mL falcon tubes (3 vials each on each
# layer) and 1mL of autoclaved groundwater was added into the vial and vortexed for 1min
# each at 1500rpm. 
# using a sterile filter we filtered the solid and got the supernatant 

# each vial was checked for ATP using 3 replicates each

# Standard processing:
conc = [50, 25, 10, 5, 1, 0.1]
measurements = [
    737494 727770 661536;
    311323 305293 338720;
    127328 124937 123190;
    54770 55178 62166;
    17081 18477 17905;
    1264 1282 1199
]
first_measurement = mean(measurements, dims=2)[:]
std_dev = std(measurements, dims=2)
std_percentage = (std_dev ./ first_measurement) * 100

# Fit linear regression
f = linear_interpolation(sort(conc), sort(first_measurement), extrapolation_bc=Line())
y_fit =f.(conc)

# Calculate R-squared
y_mean = mean(first_measurement)
SS_tot = sum((first_measurement .- y_mean).^2)
SS_res = sum((first_measurement .- y_fit).^2)
R_squared = 1 - (SS_res / SS_tot)

# Calculate ATP concentrations function
calculate_concentration = y -> f.(y)


# Processin the data from the experiment
# Top, Middle, Bottom layer data
TOP = [1146, 1050, 1064, 707, 752, 980]
MID = [2683, 2722, 2790, 3185, 3125, 3246]
BOT = [2445, 1994, 2383, 2987, 2820, 2983]

# Sample data
Col_A = [
    [1370 1406 1419 465 508 1116];
    [1241 1310 1483 1145 1087 1189];
    [5516 5546 5333 1673 1618 1589];
    [915 885 985 542 416 395]
]

Col_B = [
    [1105 1067 1091 748 1340 796];
    [743 739 832 1243 1262 1276];
    [435 459 534 503 454 454];
    [914 323 341 331 318 294]
]

Col_C = [
    [818 788 887 440 589 437];
    [1948 1942 1885 1254 1214 1214];
    [1094 1100 1007 1068 1136 1100];
    [515 572 488 506 457 469]
]

Col_D = [
    [211 208 210 191 268 204];
    [182 145 176 130 152 131];
    [88 99 90 141 150 143];
    [119 97 101 151 195 161]
]

Col_E = [
    [447 440 452 325 306 275];
    [272 598 207 159 194 163];
    [106 133 121 131 116 106];
    [278 327 216 136 137 141]
]

# Calculate average and concentrations for Column A
A_conc = calculate_concentration(Col_A)
B_conc = calculate_concentration(Col_B)
C_conc = calculate_concentration(Col_C)
D_conc = calculate_concentration(Col_D)
E_conc = calculate_concentration(Col_E)

# Solid phase DTP concentrations 
# Convert concentrations to amounts of ATP per gram of sand
volume = 1 # volume in mL
m_s = 2.5 # in g % solid mass %% **important** here in the sample it was saturated with water (hence here 2.5g is the wet mass) and at the T-zero sampling 2.5g is Dry mass (**Vitor)
# Average ATP concentration for each layer at t0
top_conc = calculate_concentration(TOP) .* volume ./ m_s
mid_conc = calculate_concentration(MID) .* volume ./ m_s
bot_conc = calculate_concentration(BOT) .* volume ./ m_s
mean_top = mean(top_conc)
neg_top = mean_top - minimum(top_conc)
pos_top = maximum(top_conc) - mean_top

mean_mid = mean(mid_conc)
neg_mid = mean_mid - minimum(mid_conc)
pos_mid = maximum(mid_conc) - mean_mid

mean_bot = mean(bot_conc)
neg_bot = mean_bot - minimum(bot_conc)
pos_bot = maximum(bot_conc) - mean_bot
top_ATP = top_conc * volume/m_s # total ATP in nM/g
mid_ATP = mid_conc * volume/m_s
bot_ATP = bot_conc * volume/m_s

# correction in the columns for the wet mass
phis = [0.279, 0.296, 0.282, 0.248, 0.262]
m_s_columns = m_s .- phis .* 1

# concentration in [nmol/g]
A_conc = A_conc ./ volume ./ m_s_columns[1]
B_conc = B_conc ./ volume ./ m_s_columns[2]
C_conc = C_conc ./ volume ./ m_s_columns[3]
D_conc = D_conc ./ volume ./ m_s_columns[4]
E_conc = E_conc ./ volume ./ m_s_columns[5]

# Calculate concentration for all data for Column A
A_mean = mean(A_conc, dims=2)
A_neg = A_mean .- minimum(A_conc, dims=2)
A_pos = maximum(A_conc, dims=2) .- A_mean
# doing the same for the rest of the columns:
B_mean = mean(B_conc, dims=2)
B_neg = B_mean .- minimum(B_conc, dims=2)
B_pos = maximum(B_conc, dims=2) .- B_mean
C_mean = mean(C_conc, dims=2)
C_neg = C_mean .- minimum(C_conc, dims=2)
C_pos = maximum(C_conc, dims=2) .- C_mean
D_mean = mean(D_conc, dims=2)
D_neg = D_mean .- minimum(D_conc, dims=2)
D_pos = maximum(D_conc, dims=2) .- D_mean
E_mean = mean(E_conc, dims=2)
E_neg = E_mean .- minimum(E_conc, dims=2)
E_pos = maximum(E_conc, dims=2) .- E_mean

# Define labels for plot
labels = ["Time zero", "Inlet", "3-6 cm", "6-9 cm", "Outlet"]
# colors for the bars
colors = [
    RGB(0.2, 0.2, 0.8); # Blue
    RGB(0.2, 0.8, 0.2); # Green
    RGB(0.8, 0.2, 0.2); # Red
    RGB(0.8, 0.8, 0.2); # Yellow
    RGB(0.6, 0.2, 0.8)  # Purple
]



# Plot each column
fig = Figure(size = (800, 600))

x = 1.5:3:12

# Sample Depth 23-24 m
ax1 = Axis(fig[1, 1], title="Sample Depth 23-24 m")
xlims!(ax1, 0, 12)
# creating the rectangle for the before experiment
pts1 = Point2f[(0, mean_mid - neg_mid), (0, mean_mid + pos_mid), (12, mean_mid + pos_mid), (12, mean_mid - neg_mid)]

er1 = errorbars!(ax1, x, vec(A_mean), vec(A_neg), vec(A_pos), label="A", color=colors[1])
sc1 = scatter!(ax1, x, vec(A_mean), color=colors[1], label="A")
l1 = lines!(ax1, x, vec(A_mean), color=colors[1], label="A", linestyle=:dash)
#poly!(Rect(0, mean_mid - neg_mid, 12, pos_mid+neg_mid), label="before experiment", color=:lightgrey, alpha=0.1, linestyle=:dash)

# Sample Depth 27-28 m
ax2 = Axis(fig[2, 1], title="Sample Depth 27-28 m", xlabel="Solid phase ATP Concentration (nM/g)")
# creating the rectangle for the before experiment
pts2 = Point2f[(0, mean_bot - neg_bot), (0, mean_bot + pos_bot), (12, mean_bot + pos_bot), (12, mean_bot - neg_bot)]

xlims!(ax2, 0, 12)
er2 = errorbars!(ax2, x, vec(B_mean), vec(B_neg), vec(B_pos), label="B", color=colors[2])
sc2 = scatter!(ax2, x, vec(B_mean), color=colors[2], label="B")
l2 = lines!(ax2, x, vec(B_mean), color=colors[2], label="B", linestyle=:dash)
er3 = errorbars!(ax2, x, vec(C_mean), vec(C_neg), vec(C_pos), label="C", color=colors[3])
sc3 = scatter!(ax2, x, vec(C_mean), color=colors[3], label="C")
l3 = lines!(ax2, x, vec(C_mean), color=colors[3], label="C", linestyle=:dash)


#oly!(ax2, pts2, label="before experiment", color=:lightgrey, alpha=0.5, linestyle=:dash)

# Link x and y axes
#linkxaxes!(ax1, ax2)
#linkyaxes!(ax1, ax2)


# Plot overall figure title
fig[:, 0] = Label(fig, "ATP Concentration Across Different Columns and Layers", fontsize=18, rotation=π/180*90)
Legend(fig[3, 1], [er1, er2, er3], ["A", "B", "C"], orientation = :horizontal, merge = true)
# Display the figure
display(fig)