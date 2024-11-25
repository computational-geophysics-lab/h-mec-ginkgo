using Plots
using DelimitedFiles

function create_plots(rootpath; skipstart=20)
    do_EVO_data = isfile(rootpath*"/EVO_data.txt")
    do_EVO_vslip = isfile(rootpath*"/EVO_vslip.txt")
    do_EVO_vxD = isfile(rootpath*"/EVO_vxD.txt")
    do_EVO_press_eff = isfile(rootpath*"/EVO_press_eff.txt")
    do_EVO_press_flu = isfile(rootpath*"/EVO_press_flu.txt")
    do_EVO_Elast_comp = isfile(rootpath*"/EVO_Elast_comp.txt")
    do_EVO_SigmaY = isfile(rootpath*"/EVO_SigmaY.txt")
    do_EVO_Sii = isfile(rootpath*"/EVO_Sii.txt")
    do_EVO_Theta = isfile(rootpath*"/EVO_Theta.txt")
    do_EVO_Visc = isfile(rootpath*"/EVO_Visc.txt")
    do_EVO_viscosity= isfile(rootpath*"/EVO_viscosity.txt")

    (!isdir(rootpath * "/plots")) && mkdir(rootpath * "/plots")

    x_fault = 1:401
    x_fault_staggered = LinRange(0.0, 400, 402)
    if do_EVO_data
        data = readdlm(rootpath * "/EVO_data.txt", skipstart = skipstart)
        p1 = plot(data[:,1]./(86400*365), data[:,3], xlabel="time [y]", ylabel="maximal slip velocity", yaxis=:log10)
        p2 = plot(data[:,3], label="slip velocity", xlabel="timestep")
        plot(p1, p2, layout=(1,2), size = (800, 400), dpi=250)
        savefig(rootpath*"/plots/EVO_data.png")
    end

    if do_EVO_vslip
        data = readdlm(rootpath * "/EVO_Vslip.txt", skipstart=skipstart)
        time = data[:,1]
        vslip = data[:, 3:end]
        heatmap(x_fault, 1:size(data)[1],data[:, 3:end], dpi=250, xlabel="Fault position [km]", ylabel="Timestep", title="Slip velocity [m/s]")
        savefig(rootpath*"/plots/EVO_vslip.png")
    end
    if do_EVO_vxD
        data = readdlm(rootpath * "/EVO_vxD.txt", skipstart=skipstart)
        time = data[:,1]
        vslip = data[:, 3:end]

        heatmap(x_fault_staggered, 1:size(data)[1],data[:, 3:end], dpi=250, xlabel="Fault position [km]", ylabel="Timestep", title="Darcy velocity vxD")
        savefig(rootpath*"/plots/EVO_vxD.png")
    end
    if do_EVO_press_eff
        data = readdlm(rootpath * "/EVO_press_eff.txt", skipstart=skipstart)
        time = data[:,1]
        vslip = data[:, 3:end]
        heatmap(x_fault_staggered, 1:size(data)[1],data[:, 3:end], dpi=250, xlabel="Fault position [km]", ylabel="Timestep", title="Effective pressure [N/m^2]")
        savefig(rootpath*"/plots/EVO_press_eff.png")
    end
    if do_EVO_press_flu
        data = readdlm(rootpath * "/EVO_press_flu.txt", skipstart=skipstart)
        time = data[:,1]
        vslip = data[:, 3:end]
        heatmap(x_fault_staggered, 1:size(data)[1],data[:, 3:end], dpi=250, xlabel="Fault position [km]", ylabel="Timestep", title="Fluid pressure [N/m^2]")
        savefig(rootpath*"/plots/EVO_press_flu.png")
    end
    if do_EVO_Elast_comp
        data = readdlm(rootpath * "/EVO_Elast_comp.txt", skipstart=skipstart)
        time = data[:,1]
        vslip = data[:, 3:end]
        heatmap(x_fault_staggered, 1:size(data)[1],data[:, 3:end], dpi=250, xlabel="Fault position [km]", ylabel="Timestep", title="Elastic compressibility")
        savefig(rootpath*"/plots/EVO_Elast_comp.png")
    end
    if do_EVO_SigmaY
        data = readdlm(rootpath * "/EVO_SigmaY.txt", skipstart=skipstart)
        time = data[:,1]
        vslip = data[:, 3:end]
        heatmap(x_fault, 1:size(data)[1],data[:, 3:end], dpi=250, xlabel="Fault position [km]", ylabel="Timestep", title="σ_y")
        savefig(rootpath*"/plots/EVO_SigmaY.png")
    end
    if do_EVO_Sii
        data = readdlm(rootpath * "/EVO_Sii.txt", skipstart=skipstart)
        time = data[:,1]
        vslip = data[:, 3:end]
        heatmap(x_fault, 1:size(data)[1],data[:, 3:end], dpi=250, xlabel="Fault position [km]", ylabel="Timestep", title="Sii")
        savefig(rootpath*"/plots/EVO_Sii.png")
    end
    if do_EVO_Theta
        data = readdlm(rootpath * "/EVO_Theta.txt", skipstart=skipstart)
        time = data[:,1]
        vslip = data[:, 3:end]
        heatmap(x_fault, 1:size(data)[1],data[:, 3:end], dpi=250, xlabel="Fault position [km]", ylabel="Timestep", title="Theta")
        savefig(rootpath*"/plots/EVO_Theta.png")
    end
    if do_EVO_Visc
        data = readdlm(rootpath * "/EVO_Visc.txt", skipstart=skipstart)
        time = data[:,1]
        vslip = data[:, 3:end]
        heatmap(x_fault_staggered, 1:size(data)[1],data[:, 3:end], dpi=250, xlabel="Fault position [km]", ylabel="Timestep", title="Visc")
        savefig(rootpath*"/plots/EVO_Visc.png")
    end
    if do_EVO_viscosity
        data = readdlm(rootpath * "/EVO_viscosity.txt", skipstart=skipstart)
        time = data[:,1]
        vslip = data[:, 3:end]
        heatmap(x_fault, 1:size(data)[1],data[:, 3:end], dpi=250, xlabel="Fault position [km]", ylabel="Timestep", title="Viscosity")
        savefig(rootpath*"/plots/EVO_viscosity.png")
    end
    return
end