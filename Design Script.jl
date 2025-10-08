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
bmax = 0.010
g = 0.008
x_p = 0.03
d = 0.015


# SIMULATION SETUP
run_simulation::Bool = true
duration = 2e-3
timestep = 5e-9
simresolution = 100
bfield_range = range(0.010,0.015,10)
mdot_range = range(0.2e-6,1e-6,10)

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

function print_performance(T, mdot, eff, Isp, Id, P, Vd)
    println("Thrust T           = $(round(T*1000, digits=3)) mN")
    println("Mass flow rate     = $(round(mdot*1e6, digits=3)) mg/s")
    println("Anode Efficiency   = $(round(eff*100, digits=1)) %")
    println("Specific Impulse   = $(round(Isp, digits=1)) s\n")
    println("Discharge Current  = $(round(Id, digits=2)) A")
    println("Input Power        = $(round(P*1e-3, digits=2)) kW")
    println("Discharge Voltage  = $(round(Vd, digits=2)) V\n")
end

function print_results(name, d, h, T, mdot, eff, Isp, r_inner, r_outer)
    println("=== $name ===")
    println("Channel diameter d = $(round(d*1000, digits=2)) mm")
    println("Channel width h    = $(round(h*1000, digits=2)) mm")
    print_performance(T, mdot, eff, Isp, P_in / V_anode, P_in, V_anode)
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

function GetMagneticField(Bmax = bmax, growth = g, decay = d, peak = x_p)
    x_bfield = 0:0.0001:(4*L_channel)
    magnetic_field = het.MagneticField(z = [], B = [])
    
    if (bfield_file != "")
        magnetic_field = het.load_magnetic_field(bfield_file)
    else
        bfield = bfield_model(Bmax, growth, decay, peak, x_bfield)
        magnetic_field = het.MagneticField(z = x_bfield, B = bfield)
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

function RunVaryMdotAndBfield(geom, simparams)
    # VARYING MASS FLOW RATE AND B FIELD SIMULATION
    sim_thrust_range = zeros(length(mdot_range), length(bfield_range))
    sim_dischargeCurrent_range = zeros(length(mdot_range), length(bfield_range))
    sim_anodeEff_range = zeros(length(mdot_range), length(bfield_range))
    sim_isp_range = zeros(length(mdot_range), length(bfield_range))
    sim_power_range = zeros(length(mdot_range), length(bfield_range))

    elements = length(mdot_range) * length(bfield_range)
    avg_time_per_sim = 20  # seconds, rough estimate
    total_estimated_time = elements * avg_time_per_sim / 60  # in minutes
    for mdot in mdot_range
        for bfield in bfield_range
            thruster = het.Thruster(name = "H-CHeT-300",
                                    geometry = geom,
                                    magnetic_field = GetMagneticField(bfield, g, d, x_p))
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
            #print_performance(sim_thrust[end], mdot, sim_anodeEff[end], sim_isp[end])
            
            sim_thrust_range[findfirst(==(mdot), mdot_range), findfirst(==(bfield), bfield_range)] = sim_thrust[end]
            sim_dischargeCurrent_range[findfirst(==(mdot), mdot_range), findfirst(==(bfield), bfield_range)] = sim_dischargeCurrent[end]
            sim_anodeEff_range[findfirst(==(mdot), mdot_range), findfirst(==(bfield), bfield_range)] = sim_anodeEff[end]
            sim_isp_range[findfirst(==(mdot), mdot_range), findfirst(==(bfield), bfield_range)] = sim_isp[end]
            sim_power_range[findfirst(==(mdot), mdot_range), findfirst(==(bfield), bfield_range)] = sim_power[end]
            #println("Completed B-field = $(bfield * 1e3) mT")
        end
        #println("Completed mdot = $(mdot * 1e6) mg/s")
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

function RunBasicParams(config, simparams)
    solution = het.run_simulation(config, simparams)

    sim_thrust, sim_dischargeCurrent, sim_anodeEff, sim_isp, sim_power = HETPerformance(solution)

    println("Final thrust: $(sim_thrust[end] * 1000) mN")
    println("Final discharge current: $(sim_dischargeCurrent[end]) A")
    println("Final anode efficiency: $(sim_anodeEff[end])")
    println("Final specific impulse: $(sim_isp[end]) s")
    println("Final input power: $(sim_power[end] * 1e-3) kW")
end

function RunVaryMdot(geom, simparams)
    # VARYING MASS FLOW RATE SIMULATION
    sim_thrust_range = zeros(length(mdot_range))
    sim_dischargeCurrent_range = zeros(length(mdot_range))
    sim_anodeEff_range = zeros(length(mdot_range))
    sim_isp_range = zeros(length(mdot_range))
    sim_power_range = zeros(length(mdot_range))

    elements = length(mdot_range)
    avg_time_per_sim = 20  # seconds, rough estimate
    total_estimated_time = elements * avg_time_per_sim / 60  # in minutes
    for mdot in mdot_range
        thruster = het.Thruster(name = "H-CHeT-300",
                                geometry = geom,
                                magnetic_field = GetMagneticField(bmax, g, d, x_p))
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
        #print_performance(sim_thrust[end], mdot, sim_anodeEff[end], sim_isp[end])
        
        sim_thrust_range[findfirst(==(mdot), mdot_range)] = sim_thrust[end]
        sim_dischargeCurrent_range[findfirst(==(mdot), mdot_range)] = sim_dischargeCurrent[end]
        sim_anodeEff_range[findfirst(==(mdot), mdot_range)] = sim_anodeEff[end]
        sim_isp_range[findfirst(==(mdot), mdot_range)] = sim_isp[end]
        sim_power_range[findfirst(==(mdot), mdot_range)] = sim_power[end]
    end
    # Plotting 1D plots for all five parameters
    fig = mk.Figure(size = (1200, 800))
    ax1 = mk.Axis(fig[1, 1], xlabel = "Mass Flow Rate (mg/s)",
                           ylabel = "Thrust (mN)",
                           title = "Thrust vs Mass Flow Rate")
    mk.lines!(ax1, mdot_range .* 1e6, sim_thrust_range .* 1000, linewidth = 2, color = :blue)
    ax2 = mk.Axis(fig[2, 1], xlabel = "Mass Flow Rate (mg/s)",
                           ylabel = "Discharge Current (A)",
                           title = "Discharge Current vs Mass Flow Rate")
    mk.lines!(ax2, mdot_range .* 1e6, sim_dischargeCurrent_range, linewidth = 2, color = :blue)
    ax3 = mk.Axis(fig[3, 1], xlabel = "Mass Flow Rate (mg/s)",
                           ylabel = "Anode Efficiency (%)",
                           title = "Anode Efficiency vs Mass Flow Rate")
    mk.lines!(ax3, mdot_range .* 1e6, sim_anodeEff_range .* 100, linewidth = 2, color = :blue)
    ax4 = mk.Axis(fig[1, 2], xlabel = "Mass Flow Rate (mg/s)",
                           ylabel = "Specific Impulse (s)",
                           title = "Specific Impulse vs Mass Flow Rate")
    mk.lines!(ax4, mdot_range .* 1e6, sim_isp_range, linewidth = 2, color = :blue)  
    ax5 = mk.Axis(fig[2, 2], xlabel = "Mass Flow Rate (mg/s)",
                           ylabel = "Input Power (kW)",
                           title = "Input Power vs Mass Flow Rate")
    mk.lines!(ax5, mdot_range .* 1e6, sim_power_range .* 1e-3, linewidth = 2, color = :blue)
    mk.display(fig)
    return
end

function RunVaryBfield(geom, simparams)
    # VARYING B FIELD SIMULATION
    sim_thrust_range = zeros(length(bfield_range))
    sim_dischargeCurrent_range = zeros(length(bfield_range))
    sim_anodeEff_range = zeros(length(bfield_range))
    sim_isp_range = zeros(length(bfield_range))
    sim_power_range = zeros(length(bfield_range))

    elements = length(bfield_range)
    avg_time_per_sim = 20  # seconds, rough estimate
    total_estimated_time = elements * avg_time_per_sim / 60  # in minutes
    for bfield in bfield_range
        thruster = het.Thruster(name = "H-CHeT-300",
                                geometry = geom,
                                magnetic_field = GetMagneticField(bfield, g, d, x_p))
        config = het.Config(thruster = thruster,
                                domain = (0.0, 4 * L_channel),
                                discharge_voltage = V_anode,
                                propellants = [het.Propellant(het.Krypton, flow_rate_kg_s = avg_mdot, max_charge = 3)])
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
        #print_performance(sim_thrust[end], mdot, sim_anodeEff[end], sim_isp[end])
        
        sim_thrust_range[findfirst(==(bfield), bfield_range)] = sim_thrust[end]
        sim_dischargeCurrent_range[findfirst(==(bfield), bfield_range)] = sim_dischargeCurrent[end]
        sim_anodeEff_range[findfirst(==(bfield), bfield_range)] = sim_anodeEff[end]
        sim_isp_range[findfirst(==(bfield), bfield_range)] = sim_isp[end]
        sim_power_range[findfirst(==(bfield), bfield_range)] = sim_power[end]
    end
    # Plotting 1D plots for all five parameters
    fig = mk.Figure(size = (1200, 800))
    ax1 = mk.Axis(fig[1, 1], xlabel = "Magnetic Field (T)",
                           ylabel = "Thrust (mN)",
                           title = "Thrust vs Magnetic Field")
    mk.lines!(ax1, bfield_range, sim_thrust_range .* 1000, linewidth = 2, color = :blue)
    ax2 = mk.Axis(fig[2, 1], xlabel = "Magnetic Field (T)",
                           ylabel = "Discharge Current (A)",
                           title = "Discharge Current vs Magnetic Field")
    mk.lines!(ax2, bfield_range, sim_dischargeCurrent_range, linewidth = 2, color = :blue)
    ax3 = mk.Axis(fig[3, 1], xlabel = "Magnetic Field (T)",
                           ylabel = "Anode Efficiency (%)",
                           title = "Anode Efficiency vs Magnetic Field")
    mk.lines!(ax3, bfield_range, sim_anodeEff_range .* 100, linewidth = 2, color = :blue)
    ax4 = mk.Axis(fig[1, 2], xlabel = "Magnetic Field (T)",
                           ylabel = "Specific Impulse (s)",
                           title = "Specific Impulse vs Magnetic Field")
    mk.lines!(ax4, bfield_range, sim_isp_range, linewidth = 2, color = :blue)  
    ax5 = mk.Axis(fig[2, 2], xlabel = "Magnetic Field (T)",
                           ylabel = "Input Power (kW)",
                           title = "Input Power vs Magnetic Field")
    mk.lines!(ax5, bfield_range, sim_power_range .* 1e-3, linewidth = 2, color = :blue)
    mk.display(fig)
    return
end

function PlotMagneticField(bfield_range, g, d, x_p)
    thruster_bfield_curves = []
    x = 0:0.0001:(4*L_channel)
    for bfield in bfield_range
        push!(thruster_bfield_curves, bfield_model(bfield, g, d, x_p, x))
    end
    return thruster_bfield_curves

    # Plot magnetic field
    fig = mk.Figure(size = (800, 600))
    ax = mk.Axis(fig[1, 1], xlabel = "Axial Position (m)",
                           ylabel = "Magnetic Field (T)",
                           title = "Magnetic Field Profile")
    for thruster_bfield in thruster_bfield_curves
        lineplot = mk.lines!(ax, x, thruster_bfield, linewidth = 2, label = "B-field $(thruster_bfield.B[findmax(thruster_bfield.B)[2]]*1e3) mT")
    end
    vlineplot = mk.vlines!(ax, L_channel, color = :red, label = "Channel Start/End")    
    mk.Legend(fig[1, 2], [lineplot, vlineplot], ["B-field", "Channel Start/End"], orientation = :vertical)
    mk.display(fig)
    return
end



function PlotPlasmaProperties(solution, name)
    avg = het.time_average(solution)
    z_cm = avg.grid.cell_centers .* 100

    if name == "Argon"
        # neutrals
        nn = avg.frames[].neutrals[:Ar].n
        un = avg.frames[].neutrals[:Ar].u ./ 1000 # km/s

        # ions
        ni = avg.frames[].ions[:Ar][1].n
        ui = avg.frames[].ions[:Ar][1].u ./ 1000 # km/s
        ni2 = avg.frames[].ions[:Ar][2].n
        ui2 = avg.frames[].ions[:Ar][2].u ./ 1000
        ni3 = avg.frames[].ions[:Ar][3].n
        ui3 = avg.frames[].ions[:Ar][3].u ./ 1000
    elseif name == "Xenon"
        # neutrals
        nn = avg.frames[].neutrals[:Xe].n
        un = avg.frames[].neutrals[:Xe].u ./ 1000 # km/s

        # ions
        ni = avg.frames[].ions[:Xe][1].n
        ui = avg.frames[].ions[:Xe][1].u ./ 1000 # km/s
        ni2 = avg.frames[].ions[:Xe][2].n
        ui2 = avg.frames[].ions[:Xe][2].u ./ 1000
        ni3 = avg.frames[].ions[:Xe][3].n
        ui3 = avg.frames[].ions[:Xe][3].u ./ 1000
    elseif name == "Krypton"
        # neutrals
        nn = avg.frames[].neutrals[:Kr].n
        un = avg.frames[].neutrals[:Kr].u ./ 1000 # km/s

        # ions
        ni = avg.frames[].ions[:Kr][1].n
        ui = avg.frames[].ions[:Kr][1].u ./ 1000 # km/s
        ni2 = avg.frames[].ions[:Kr][2].n
        ui2 = avg.frames[].ions[:Kr][2].u ./ 1000
        ni3 = avg.frames[].ions[:Kr][3].n
        ui3 = avg.frames[].ions[:Kr][3].u ./ 1000
    else
        error("Unsupported propellant: $name")
    end

    # electrons
    ne = avg.frames[].ne
    Te = avg.frames[].Tev
    ue = avg.frames[].ue ./ 1000 # electron velocity in km/s

    # Field parameters normalised from 0 to 1
    B = avg.frames[].B
    psi = avg.frames[].potential
    B ./= maximum(B)
    psi ./= maximum(psi)

    f = mk.Figure(size = (1000, 700))

    # 2x2 layout
    ax1 = mk.Axis(f[1, 1], xlabel = "Axial position [cm]", ylabel = "Density [m^-3]", title = "Channel species densities ($name)", yscale = mk.log10)
    l1 = mk.lines!(ax1, z_cm, nn, label = "Neutrals")
    l3 = mk.lines!(ax1, z_cm, ne, label = "Electrons")
    l2 = mk.lines!(ax1, z_cm, ni, label = "Singly charged ions")
    l4 = mk.lines!(ax1, z_cm, ni2, label = "Doubly charged ions")
    l5 = mk.lines!(ax1, z_cm, ni3, label = "Triply charged ions")
    mk.Legend(f[1, 2], [l1, l2, l3, l4, l5], ["Neutrals", "Singly charged ions", "Electrons", "Doubly charged ions", "Triply charged ions"]; orientation = :vertical)

    ax2 = mk.Axis(f[2, 1], xlabel = "Axial position [cm]", ylabel = "Velocity [km/s]", title = "Channel species velocities ($name)")
    l6 = mk.lines!(ax2, z_cm, un, label = "Neutrals")
    l10 = mk.lines!(ax2, z_cm, ue, label = "Electrons") 
    l7 = mk.lines!(ax2, z_cm, ui, label = "Singly charged ions")
    l8 = mk.lines!(ax2, z_cm, ui2, label = "Doubly charged ions")
    l9 = mk.lines!(ax2, z_cm, ui3, label = "Triply charged ions")
    mk.Legend(f[2, 2], [l6, l7, l8, l9, l10], ["Neutrals", "Singly charged ions", "Doubly charged ions", "Triply charged ions", "Electrons"]; orientation = :vertical)

    ax3 = mk.Axis(f[1, 3], xlabel = "Axial position [cm]", ylabel = "Energy [eV]", title = "Electron temperature ($name)")
    mk.lines!(ax3, z_cm, Te)

    ax4 = mk.Axis(f[2, 3], xlabel = "Axial position [cm]", ylabel = "Normalised Field", title = "Electric and Magnetic fields ($name)")
    l12 = mk.lines!(ax4, z_cm, B, label = "Magnetic field")
    l13 = mk.lines!(ax4, z_cm, psi, label = "Potential")
    mk.Legend(f[2, 4], [l12, l13], ["Magnetic field", "Potential"]; orientation = :vertical)

    mk.display(f)
end

# ADD A CUSTOM INPUTS SIMULATION AS WELL AS A SCALING LAW SIMULATION
if run_simulation == true
    println("Starting simulation...")

    # SINGLE POINT SIMULATION
    # === SIMULATION ===
    geom = het.Geometry1D(channel_length = L_channel,
                            inner_radius = avg_r_inner,
                            outer_radius = avg_r_outer)
    thruster = het.Thruster(name = "H-CHeT-300",
                            geometry = geom,
                            magnetic_field = GetMagneticField(bmax, g, d, x_p))
    config = het.Config(thruster = thruster,
                            domain = (0.0, 3 * L_channel),
                            discharge_voltage = V_anode,
                            propellants = [het.Propellant(het.Krypton, flow_rate_kg_s = avg_mdot, max_charge = 3)])
    simparams = het.SimParams(grid = het.UnevenGrid(simresolution),
                            dt = timestep,
                            duration = duration,
                            num_save = Int(duration * 1e6))
    #PlotMagneticField(bmax, g, d, x_p)
    RunVaryMdot(geom, simparams)
    RunVaryBfield(geom, simparams)
    #RunBasicParams(config, simparams)
    #RunVaryMdotAndBfield(geom, simparams)
    solution = het.run_simulation(config, simparams)
    PlotPlasmaProperties(solution, "Krypton")
    thrust, dischargeCurrent, anodeEff, isp, power = HETPerformance(solution)
    
    print_performance(thrust[end], avg_mdot[end], anodeEff[end], isp[end], dischargeCurrent[end], power[end], V_anode)
    return
end
