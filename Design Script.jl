using HallThruster: HallThruster as het
using CairoMakie: Makie as mk

# CONSTANTS
# PHYSICS
m_e = 9.1093837e-31
k_b = 1.380649e-23
# LEE
C_T_Lee = 892.7
C_P_Lee = 633
C_mdot_Lee = 0.003
C_hd_Lee = 0.242

# Dannenmeyer
C_T1_Dann  = 1077.3
C_T2_Dann  = 0.924
C_Isp_Dann = 109.835
C_P_Dann   = 11928
C_hd_Dann  = 0.242

# Hybrid
C_T1_Hyb   = 854.1
C_T2_Hyb   = 0.7153
C_Isp_Hyb  = 86.3
C_P_Hyb    = 11711
C_hd_Hyb   = 0.256

# INPUTS
P_in = 100
V_anode = 300
L_channel = 0.03
bfield_file::String = ""
bmax = 0.02
g = 0.01
x_p = 0.03
d = 0.04


# SIMULATION SETUP
run_simulation::Bool = true
duration = 2e-3
timestep = 1e-9
simresolution = 100
bfield_range = range(0.01,0.02,10)
mdot_range = range(0.4e-6,1e-6,10)

# === SCALING MATHS ===
# DANNENMEYER MODEL
function dannenmeyer_model(P, Ud)
    d = sqrt(P / (sqrt(Ud) * C_P_Dann))
    h = d * C_hd_Dann
    T = C_T2_Dann * d^2 * sqrt(Ud)
    mdot = T / (C_T1_Dann * sqrt(Ud))
    return (d, h, T, mdot)
end

# LEE MODEL
function lee_model(P, Ud)
    d = sqrt(P / (C_P_Lee * Ud))
    h = C_hd_Lee * d
    mdot = C_mdot_Lee * h * d
    T = C_T_Lee * mdot * sqrt(Ud)
    return (d, h, T, mdot)
end

# HYBRID MODEL (dannenmeyer equations using lee database)
function hybrid_model(P, Ud)
    d = sqrt(P / (sqrt(Ud) * C_P_Hyb))
    h = d * C_hd_Hyb
    T = C_T2_Hyb * d^2 * sqrt(Ud)
    mdot = T / (C_T1_Hyb * sqrt(Ud))
    return (d, h, T, mdot)
end

DannOutputs = dannenmeyer_model(P_in, V_anode)
LeeOutputs = lee_model(P_in, V_anode)
HybridOutputs = hybrid_model(P_in, V_anode)

function calculate_anode_efficiency(T, P, mdot)
    return T^2 / (2 * mdot * P)
end

function calculate_isp(eff, P, mdot)
    Pjet = eff * P
    vex = sqrt(Pjet * 2 / mdot)
    return vex / 9.81
end

function calculate_inner_diameter(d, h)
    return -h + sqrt(d^2 - h^2)
end

function calculate_outer_diameter(d, h)
    return h + sqrt(d^2 - h^2)
end

# OUTPUTS
Dann_d, Dann_h, Dann_T, Dann_mdot = dannenmeyer_model(P_in, V_anode)
Lee_d, Lee_h, Lee_T, Lee_mdot = lee_model(P_in, V_anode)
Hyb_d, Hyb_h, Hyb_T, Hyb_mdot = hybrid_model(P_in, V_anode)
avg_d = (Dann_d + Lee_d + Hyb_d) / 3
avg_h = (Dann_h + Lee_h + Hyb_h) / 3
avg_T = (Dann_T + Lee_T + Hyb_T) / 3
avg_mdot = (Dann_mdot + Lee_mdot + Hyb_mdot) / 3
# Calculate efficiencies
Dann_eff = calculate_anode_efficiency(Dann_T, P_in, Dann_mdot)
Lee_eff  = calculate_anode_efficiency(Lee_T,  P_in, Lee_mdot)
Hyb_eff  = calculate_anode_efficiency(Hyb_T,  P_in, Hyb_mdot)
avg_eff = (Dann_eff + Lee_eff + Hyb_eff) / 3

# Calculate inner and outer channel radius
Dann_r_inner = calculate_inner_diameter(Dann_d, Dann_h) / 2
Dann_r_outer = calculate_outer_diameter(Dann_d, Dann_h) / 2
Lee_r_inner  = calculate_inner_diameter(Lee_d, Lee_h) / 2
Lee_r_outer  = calculate_outer_diameter(Lee_d, Lee_h) / 2
Hyb_r_inner  = calculate_inner_diameter(Hyb_d, Hyb_h) / 2
Hyb_r_outer  = calculate_outer_diameter(Hyb_d, Hyb_h) / 2
avg_r_inner = (Dann_r_inner + Lee_r_inner + Hyb_r_inner) / 3
avg_r_outer = (Dann_r_outer + Lee_r_outer + Hyb_r_outer) / 3

# Calculate specific impulses
Dann_isp = calculate_isp(Dann_eff, P_in, Dann_mdot)
Lee_isp  = calculate_isp(Lee_eff,  P_in, Lee_mdot)
Hyb_isp  = calculate_isp(Hyb_eff,  P_in, Hyb_mdot)
avg_isp = (Dann_isp + Lee_isp + Hyb_isp) / 3

function print_performance(T, mdot, eff, Isp)
    println("Thrust T           = $(round(T*1000, digits=3)) mN")
    println("Mass flow rate     = $(round(mdot*1e6, digits=3)) mg/s")
    println("Anode Efficiency   = $(round(eff*100, digits=1)) %")
    println("Specific Impulse   = $(round(Isp, digits=1)) s\n")
end

function print_results(name, d, h, T, mdot, eff, Isp, r_inner, r_outer)
    println("=== $name ===")
    println("Channel diameter d = $(round(d*1000, digits=2)) mm")
    println("Channel width h    = $(round(h*1000, digits=2)) mm")
    print_performance(T, mdot, eff, Isp)
    println("Inner channel radius = $(round(r_inner*1000, digits=2)) mm")
    println("Outer channel radius = $(round(r_outer*1000, digits=2)) mm\n")
end



# Print results for all three models
print_results("Dannenmeyer & Mazouffre", Dann_d, Dann_h, Dann_T, Dann_mdot, Dann_eff, Dann_isp, Dann_r_inner, Dann_r_outer)
print_results("Lee et al.", Lee_d, Lee_h, Lee_T, Lee_mdot, Lee_eff, Lee_isp, Lee_r_inner, Lee_r_outer)
print_results("Hybrid", Hyb_d, Hyb_h, Hyb_T, Hyb_mdot, Hyb_eff, Hyb_isp, Hyb_r_inner, Hyb_r_outer)
print_results("Averaged", avg_d, avg_h, avg_T, avg_mdot, avg_eff, avg_isp, avg_r_inner, avg_r_outer)


function bfield_model(B_max, g, d, x_p, x)
B = similar(x)  # allocate result with same shape as x

    for (i, xi) in pairs(x)
        if xi >= 0 && xi < x_p
            B[i] = B_max * exp(-((xi - x_p)^2) / (2 * g^2))
        elseif xi >= x_p
            B[i] = B_max * exp(-((xi - x_p)^2) / (2 * d^2))
        end
    end
    return B
end

function GetMagnetiicField(Bmax = bmax, growth = g, decay = d, peak = x_p)
    x_bfield = 0:0.0001:(4*L_channel)
    magnetic_field = het.MagneticField(z = [], B = [])
    
    if (bfield_file != "")
        magnetic_field = het.load_magnetic_field(bfield_file)
    else
        bfield = bfield_model(Bmax, growth, decay, peak, x_bfield)
        magnetic_field = het.MagneticField(z = x_bfield, B = bfield)
    
        # # Plot magnetic field
        # fig = mk.Figure(size = (800, 600))
        # ax = mk.Axis(fig[1, 1], xlabel = "Axial Position (m)",
        #                        ylabel = "Magnetic Field (T)",
        #                        title = "Magnetic Field Profile")
        # lineplot = mk.lines!(ax, x_bfield, b_field, color = :blue, linewidth = 2, label = "B-field")
        # vlineplot = mk.vlines!(ax, L_channel, color = :red, label = "Channel Start/End")
        # mk.Legend(fig[1, 2], [lineplot, vlineplot], ["B-field", "Channel Start/End"], orientation = :vertical)
    
        # mk.display(fig)
    end
    return magnetic_field
end

function HETPerformance(solution)
    thrust = het.thrust(solution)
    dischargeCurrent = het.discharge_current(solution)
    anodeEff = het.anode_eff(solution)
    isp = thrust ./ (avg_mdot * 9.80665)
    power = V_anode .* dischargeCurrent

    return thrust, dischargeCurrent, anodeEff, isp, power
end

if run_simulation == true
    println("Starting simulation...")

    # SINGLE POINT SIMULATION
    # === SIMULATION ===
    geom = het.Geometry1D(channel_length = L_channel,
                            inner_radius = avg_r_inner,
                            outer_radius = avg_r_outer)
    thruster = het.Thruster(name = "H-CHeT-300",
                            geometry = geom,
                            magnetic_field = GetMagnetiicField(bmax, g, d, x_p))
    config = het.Config(thruster = thruster,
                            domain = (0.0, 4 * L_channel),
                            discharge_voltage = V_anode,
                            propellants = [het.Propellant(het.Krypton, flow_rate_kg_s = avg_mdot, max_charge = 3)])
    simparams = het.SimParams(grid = het.EvenGrid(simresolution),
                            dt = timestep,
                            duration = duration,
                            num_save = Int(duration * 1e6))
    solution = het.run_simulation(config, simparams)

    sim_thrust, sim_dischargeCurrent, sim_anodeEff, sim_isp, sim_power = HETPerformance(solution)

    println("Final thrust: $(sim_thrust[end] * 1000) mN")
    println("Final discharge current: $(sim_dischargeCurrent[end]) A")
    println("Final anode efficiency: $(sim_anodeEff[end])")
    println("Final specific impulse: $(sim_isp[end]) s")
    println("Final input power: $(sim_power[end] * 1e-3) kW")

    # VARYING MASS FLOW RATE AND B FIELD SIMULATION
    sim_thrust_range = zeros(length(mdot_range), length(bfield_range))
    sim_dischargeCurrent_range = zeros(length(mdot_range), length(bfield_range))
    sim_anodeEff_range = zeros(length(mdot_range), length(bfield_range))
    sim_isp_range = zeros(length(mdot_range), length(bfield_range))
    sim_power_range = zeros(length(mdot_range), length(bfield_range))

    for mdot in mdot_range
        for bfield in bfield_range
            thruster = het.Thruster(name = "H-CHeT-300",
                                    geometry = geom,
                                    magnetic_field = GetMagnetiicField(bfield, g, d, x_p))
            config = het.Config(thruster = thruster,
                                    domain = (0.0, 4 * L_channel),
                                    discharge_voltage = V_anode,
                                    propellants = [het.Propellant(het.Krypton, flow_rate_kg_s = mdot, max_charge = 3)])
            solution = het.run_simulation(config, simparams)

            sim_thrust, sim_dischargeCurrent, sim_anodeEff, sim_isp, sim_power = HETPerformance(solution)

            if (sim_thrust[end] < 0) || (sim_dischargeCurrent[end] < 0) || (sim_anodeEff[end] < 0) || (sim_isp[end] < 0) || (sim_power[end] < 0)
                sim_thrust[end] = 0
                sim_dischargeCurrent[end] = 0
                sim_anodeEff[end] = 0
                sim_isp[end] = 0
                sim_power[end] = 0
                continue
            end
            print_performance(sim_thrust[end], mdot, sim_anodeEff[end], sim_isp[end])
            
            sim_thrust_range[findfirst(==(mdot), mdot_range), findfirst(==(bfield), bfield_range)] = sim_thrust[end]
            sim_dischargeCurrent_range[findfirst(==(mdot), mdot_range), findfirst(==(bfield), bfield_range)] = sim_dischargeCurrent[end]
            sim_anodeEff_range[findfirst(==(mdot), mdot_range), findfirst(==(bfield), bfield_range)] = sim_anodeEff[end]
            sim_isp_range[findfirst(==(mdot), mdot_range), findfirst(==(bfield), bfield_range)] = sim_isp[end]
            sim_power_range[findfirst(==(mdot), mdot_range), findfirst(==(bfield), bfield_range)] = sim_power[end]
            println("Completed B-field = $(bfield * 1e3) mT")
        end
        println("Completed mdot = $(mdot * 1e6) mg/s")
    end

    # Plotting contours for all five parameters
    fig = mk.Figure(size = (1200, 800))
    ax1 = mk.Axis(fig[1, 1], xlabel = "Magnetic Field (T)",
                           ylabel = "Mass Flow Rate (mg/s)",
                           title = "Thrust (mN)")
    co = mk.contourf!(ax1, bfield_range, mdot_range .* 1e6, sim_thrust_range .* 1000, colormap = :viridis)
    mk.Colorbar(fig[1, 2], co, label = "Thrust (mN)")

    ax2 = mk.Axis(fig[2, 1], xlabel = "Magnetic Field (T)",
                           ylabel = "Mass Flow Rate (mg/s)",
                           title = "Discharge Current (A)")
    co2 = mk.contourf!(ax2, bfield_range, mdot_range .* 1e6, sim_dischargeCurrent_range, colormap = :viridis)
    mk.Colorbar(fig[2, 2], co2, label = "Discharge Current (A)")

    ax3 = mk.Axis(fig[1, 3], xlabel = "Magnetic Field (T)",
                           ylabel = "Mass Flow Rate (mg/s)",
                           title = "Anode Efficiency")

    co3 = mk.contourf!(ax3, bfield_range, mdot_range .* 1e6, sim_anodeEff_range .* 100, colormap = :viridis)
    mk.Colorbar(fig[1, 4], co3, label = "Anode Efficiency (%)")
    
    ax4 = mk.Axis(fig[2, 3], xlabel = "Magnetic Field (T)",
                           ylabel = "Mass Flow Rate (mg/s)",
                           title = "Specific Impulse (s)")
    co4 = mk.contourf!(ax4, bfield_range, mdot_range .* 1e6, sim_isp_range, colormap = :viridis)
    mk.Colorbar(fig[2, 4], co4, label = "Specific Impulse (s)")

    ax5 = mk.Axis(fig[3, 1], xlabel = "Magnetic Field (T)",
                           ylabel = "Mass Flow Rate (mg/s)",
                           title = "Input Power (kW)")

    co5 = mk.contourf!(ax5, bfield_range, mdot_range .* 1e6, sim_power_range .* 1e-3, colormap = :viridis)
    mk.Colorbar(fig[3, 2], co5, label = "Input Power (kW)")

    mk.display(fig)
    return
end
