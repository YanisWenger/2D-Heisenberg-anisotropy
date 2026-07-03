using Distributions, Plots, Makie, CairoMakie, ColorSchemes, Colors, CSV, DataFrames, Dierckx, LsqFit, JLD2, Glob, MCMCChains   # To have basic math (mean), To plot, To plot, To plot (needed for Makie), colors for lattices, write/read csv (CSV & DataFrames), interpolate & fit, save data in compact file, glob to read files and name files, To check convergence
include("./Heisenberg_analysis_functions.jl")


# ----- Initial parameters ----- #


const N    = 1000_000     # Number of lattice sweeps
const L    = [300]  # Lattices size
const Nperbin = 100       # Number of measurements per bins
const PBC  = true         # Periodic Boundary Conditions
folder_name_end=""    # if the folder as a name like [8, 12, 20]_1000000_Ising, the name end is _Ising. By default the folders are created without name end


# ----- Get data ----- #

# all quantities are intensive (E and M are saved normalized, C and χ are then computed directly normalized)

const T, D_for_L, T_for_LD = Get_T_LD(N, L, folder_name_end)   # Obtain Temperature and anisotropic term in Data/$(L)_$N i.e. Data/[8, 12, 20]_1000000
Nmeasurement = length(load("Data/$(L)_$N$(folder_name_end)/$(L[1])_$(T_for_LD[(L[1],D_for_L[L[1]][1])][1])_$(D_for_L[L[1]][1]).jld2")[:"Energies"])
Nbin = div(Nmeasurement, Nperbin)

nL, nT, nD  =   length(L), length(T), maximum([length(D_for_L[i]) for i in L]) # number of different lattice sizes, temperatures and anisotropic term
E, EΔ       =   zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD)
M, MΔ       =   zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD)
Mx, My, Mz, MxΔ, MyΔ, MzΔ = zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD)
χ, χΔ       =   zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD)
χx, χy, χz, χxΔ, χyΔ, χzΔ       =   zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD)
C, CΔ #=, CΔ2=#  =   zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD)#, zeros(Float32, nL, nT, nD)
Acceptance  =   zeros(4, nL, nT, nD)
Corr        =   Array{Vector{Float32}}(undef, nL, nT, nD)
# AllLattices =   Array{Vector{Array{Float32, 3}}}(undef, nL, nT, nD) # for Correlation time (do not exists by default)
# SwapAccept  =   Array{Vector}(undef, nL, nD)    # acceptance of the lattices swap (parallel tempering)

for l in eachindex(L)
    ll = L[l]
    for d in eachindex(D_for_L[ll])
        dd = D_for_L[ll][d]
        for t in eachindex(T_for_LD[(ll, dd)])
            tt = T_for_LD[(ll, dd)][t]
            Data = load("Data/$(L)_$N$(folder_name_end)/$(ll)_$(tt)_$(dd).jld2")
            En = Data[:"Energies"]
            mx = Data[:"Mx"]
            my = Data[:"My"]
            mz = Data[:"Mz"]
            Ma = sqrt.(mx.^2+my.^2+mz.^2)
            E[l,t,d], EΔ[l,t,d] = mean(En),                             std(Binor(En, Nbin, Nperbin))/sqrt(Nbin)
            M[l,t,d], MΔ[l,t,d] = mean(Ma),                             std(Binor(Ma, Nbin, Nperbin))/sqrt(Nbin)
            Mx[l,t,d], MxΔ[l,t,d] = mean(mx),                           std(Binor(mx, Nbin, Nperbin))/sqrt(Nbin)
            My[l,t,d], MyΔ[l,t,d] = mean(my),                           std(Binor(my, Nbin, Nperbin))/sqrt(Nbin)
            Mz[l,t,d], MzΔ[l,t,d] = mean(mz),                           std(Binor(mz, Nbin, Nperbin))/sqrt(Nbin)
            χ[l,t,d], χΔ[l,t,d] = (mean(Ma.^2)-M[l,t,d]^2)/tt*ll^2, std(Binor2nd(Ma, Nbin, Nperbin, tt, ll))/sqrt(Nbin)
            χx[l,t,d], χxΔ[l,t,d] = (mean(mx.^2)-Mx[l,t,d]^2)/tt*ll^2, std(Binor2nd(mx, Nbin, Nperbin, tt, ll))/sqrt(Nbin)
            χy[l,t,d], χyΔ[l,t,d] = (mean(my.^2)-My[l,t,d]^2)/tt*ll^2, std(Binor2nd(my, Nbin, Nperbin, tt, ll))/sqrt(Nbin)
            χz[l,t,d], χzΔ[l,t,d] = (mean(mz.^2)-Mz[l,t,d]^2)/tt*ll^2, std(Binor2nd(mz, Nbin, Nperbin, tt, ll))/sqrt(Nbin)
            C[l,t,d], CΔ[l,t,d]#=, CΔ2[l,t,d]=# = (mean(En.^2)-E[l,t,d]^2)/tt^2*ll^2, std(Binor2nd(En, Nbin, Nperbin, tt, ll))/tt/sqrt(Nbin)#, Errorpropagation(Binor(En, Nbin, Nperbin), EΔ[l,t,d])/tt^2*ll^2
            Corr[l,t,d]         = Data[:"corr"]
            Acceptance[:,l,t,d] = Data[:"accept"]
            # AllLattices[l,t,d] = Data[:"Lattices"]
            # println("N = ", ll, "\tT = ", tt, " \td = ", dd, "\t E = ", round(E[l,t,d];digits=3), " ± ", round(EΔ[l,t,d];digits=5), "\tM = ", round(M[l,t,d];digits=3), " ± ", round(MΔ[l,t,d];digits=5), " \t χ = ", round(χ[l,t,d];digits=3), " ± ", round(χΔ[l,t,d];digits=3), "\tC = ", round(C[l,t,d];digits=3), " ± ", round(CΔ[l,t,d];digits=4), "\taccept = ", round.(Acceptance[:,l,t,d];digits=3))
        end
        # SwapAccept[l,d] = load("Data/$(L)_$N$(folder_name_end)/swap_$(ll)_$(dd)_$(T_for_LD[(ll, dd)][1]).jld2")[:"SwapAccept"]
    end
    print(ll, "     ")
end


# X = load("AllData.jld2")
# E=X["E"]; EΔ=X["EΔ"]; C=X["C"]; CΔ=X["CΔ"]; M=X["M"]; MΔ=X["MΔ"]; Mx=X["Mx"]; MxΔ=X["MxΔ"]; My=X["My"]; MyΔ=X["MyΔ"]; Mz=X["Mz"]; MzΔ=X["MzΔ"]; χ=X["χ"]; χΔ=X["χΔ"]; χx=X["χx"]; χxΔ=X["χxΔ"]; χy=X["χy"]; χyΔ=X["χyΔ"]; χz=X["χz"]; χzΔ=X["χzΔ"]; Corr=X["Corr"]; const T_for_LD=X["T_for_LD"]; const L=X["L"]; const D_for_L=X["D_for_L"]



# -----     Plots and analysis     ----- #


PlotColors = cgrad([RGB(.8,.9,1), RGB(0,0,.5), RGB(0,0,0), RGB(.6,0,0), RGB(1,.75,.75)], [0., .4999, .5, .5001, 1.], categorical = false) # Color for plots
colorbar(PlotColors, D) # Need to create the colorbar for some following graphs    

plotL(L, T_for_LD, D_for_L, M,- 1, "M", MΔ)  # Plot all lattice sizes for the same d.    Fifth entry for which d
plotL(L, T_for_LD, D_for_L, Mx, -1, "M_x", MxΔ)
plotL(L, T_for_LD, D_for_L, χ, .7, "Susc", χΔ)
plotL(L, T_for_LD, D_for_L, χx, -1, "Susc_x", χxΔ)
plotL(L, T_for_LD, D_for_L, χy, -1, "Susc_y", χyΔ)
plotL(L, T_for_LD, D_for_L, χz, -1, "Susc_z", χzΔ)
plotL(L, T_for_LD, D_for_L, E, 1, "E", EΔ)
plotL(L, T_for_LD, D_for_L, C, .7, "C", CΔ)

plotD(PlotColors, L, T_for_LD, D_for_L, χ, "Susc", χΔ)   # Plot all the d for a given lattice size (last number for lattice size : L[#], the bigger one by default)
# plotD(L, T_for_LD, D, C, "C", CΔ2, 3)   # with error propag from ΔE, smaller error bars
plotD(PlotColors, L, T_for_LD, D_for_L, C, "C", CΔ)
plotD(PlotColors, L, T_for_LD, D_for_L, E, "E", EΔ, 1)
plotD(PlotColors, L, T_for_LD, D_for_L, M, "M", MΔ)

plotTCorrelation(L, T_for_LD, D_for_L, Corr, 1, 4)
plotLCorrelation(L, T_for_LD, D_for_L, Corr, 1, .8)
plotDCorrelation(L, T_for_LD, D_for_L, Corr, 4, .8)


# crit(L, T_for_LD[(L[1], D[1])], C[:,:,1])
# α, σα, γ, σγ = zeros(nD), zeros(nD), zeros(nD), zeros(nD)
# for i=1:nD
#     α[i], σα[i] = crit(L, T, C[:,:,i])#, "Capacity")
#     γ[i], σγ[i] = crit(L, T, χ[:,:,i])#, "Susceptibility")
# end
# println(α, "\n", σα, "\n\n", γ, "\n", σγ)



ξ, σξ, p = critlength(Corr, T_for_LD, L, D_for_L, 9, 0.7, 1.0, true, false)           # sixth entry for the wanted temperature, seventh to plot or not, eighth component to "true" to get exponent of the algebraic decay instead of correlation length

let
l=9; d=-.5
ξ, σξ, p1 = critlength(Corr, T_for_LD, L, D_for_L, l, d, 0.0, false, false)           # sixth entry for the wanted temperature, seventh to plot or not, eighth component to "true" to get exponent of the algebraic decay instead of correlation length
ξ, σξ, p2 = critlength(Corr, T_for_LD, L, D_for_L, l, d, 0.0, false, true)           # sixth entry for the wanted temperature, seventh to plot or not, eighth component to "true" to get exponent of the algebraic decay instead of correlation length
ξ, σξ, p3 = critlength(Corr, T_for_LD, L, D_for_L, l, d, 10, false, false)           # sixth entry for the wanted temperature, seventh to plot or not, eighth component to "true" to get exponent of the algebraic decay instead of correlation length
ξ, σξ, p4 = critlength(Corr, T_for_LD, L, D_for_L, l, d, 10, false, true)           # sixth entry for the wanted temperature, seventh to plot or not, eighth component to "true" to get exponent of the algebraic decay instead of correlation length
p=Plots.plot(p1, p2, p3, p4, layout=(2, 2), size=(800, 600), sgtitle="Correlation for L=$(L[l]), D=$(d)")
display(p)
# Plots.savefig("Plot/Correlation_$(L[l])_$(d).pdf")
end



Plot_Max_ξ(L, D, T_for_LD, Corr, 7, true) # bad for now as it always tries to fit an exponential, but below Tc, this is a power law

# CorrTimePlot(AllLattices, T_for_LD, D, 3, L, 2, false)   # Works only when lattices have been saved in the simulation process, doesn't by default.   Fourth element for which d, sixth for which L

histogram([x for cell in SwapAccept for x in cell])


Tmax_χ, χmax, χfitln, χfitln_err, χfitpower, χfit_Tc = Plot_Max_C_Χ(PlotColors, D_for_L, T_for_LD, χ, χΔ, L, "Susc") # Plot and fit the max of C or χ
Tmax_C, Cmax, Cfitln, Cfitln_err, Cfitpower, Cfit_Tc = Plot_Max_C_Χ(PlotColors, D_for_L, T_for_LD, C, CΔ, L, "C")


Plots.plot(D, [x[1][1] for x in χfitpower], seriestype=:scatter, yerr=[x[2][1] for x in χfitpower], label="", title="χ ∝ L^-γ vs D", xlabel="⬅ Ising                              Value of D                              XY ➡", ylabel="-γ")
Plots.plot(D, [x[1][1] for x in Cfitpower], seriestype=:scatter, yerr=[x[2][1] for x in Cfitpower], label="", title="C ∝ L^a vs D", xlabel="⬅ Ising                              Value of D                              XY ➡", ylabel="a")

Plots.plot(D, [x[1] for x in χfitln], seriestype=:scatter, yerr=[x[2] for x in χfitln], label="", title="χ = a*ln(L) vs D", xlabel="⬅ Ising                              Value of D                              XY ➡", ylabel="a")
Plots.plot(D, Cfitln[:,1], seriestype=:scatter, yerr=Cfitln_err[:,1], label="", title="C = a*ln(L/L0)*(1+L'/L) vs D", xlabel="⬅ Ising                              Value of D                              XY ➡", ylabel="a")
Plots.plot(D, Cfitln[:,2], seriestype=:scatter, yerr=Cfitln_err[:,2], label="", title="C = a*ln(L/L0)*(1+L'/L) vs D", xlabel="⬅ Ising                              Value of D                              XY ➡", ylabel="L0")
Plots.plot(D, Cfitln[:,3], seriestype=:scatter, yerr=Cfitln_err[:,3], label="", title="C = a*ln(L/L0)*(1+L'/L) vs D", xlabel="⬅ Ising                              Value of D                              XY ➡", ylabel="L'")

Plots.plot(D, [x[1][2] for x in χfit_Tc], seriestype=:scatter, yerr=[x[2][2] for x in χfit_Tc], label="", title = "T of χ_max for L → ∞ vs D", xlabel="⬅ Ising                              Value of D                              XY ➡", ylabel="T of χ_max")


scaling_plot_C_Ising(L, T_for_LD, D_for_L, C, 0.7, Tmax_C, CΔ)
scaling_plot_χ_Ising(L, T_for_LD, D_for_L, χ, 0.7, Tmax_χ, χΔ)
Plots.savefig("Plot/C_aln(L).pdf")


# fit corrlength (see where critical length diverges)
# correlation vs distance is not correlation length
# Yeomans 116 (115-118) to rescale C and 


for d=1:23
    p=Plots.plot(T_for_LD[L[1],D[d]], χ[1,1:length(T_for_LD[L[1],D[d]]),d], label="χ $(D[d])")
    display(p)
    p=Plots.plot(T_for_LD[L[1],D[d]], C[1,1:length(T_for_LD[L[1],D[d]]),d], label="C $(D[d])")
    display(p)
end

T_for_LD[L[1],D[20]][10:20]


Tpic=Vector{Vector{Float32}}(undef, nD)
l=1
v = vcat(collect(-.004f0:.001f0:-.001f0), collect(.001f0:.001f0:.004f0))
for d=1:23
    # Tpic[d-10]=Float32.(round.(T_for_LD[(L[l], D[d])][argmax(χ[l,:,d])] .+ v;digits=3))
    Tpic[d]=Float32.(round.(unique(vcat(T_for_LD[(L[l], D[d])][argmax(C[l,:,d])] .+ v,    T_for_LD[(L[l], D[d])][argmax(χ[l,:,d])] .+ v));digits=3))
end
@save "Tpic_$(L[l]).jld2" Tpic D=D[1:end]

Plots.plot(T_for_LD[L[1],D[1]], C[1,1:length(T_for_LD[L[1],D[1]]),1], label="C")
T_for_LD[L[1],D[1]]
let
    d=12
    println(D[d])
    for l=4:8
        println("$(L[l]) \t $(Tmax_C[d,l]) \t $(Tmax_χ[d,l])")
    end
    plotL(L, T_for_LD, D_for_L, C, D[d], "C", CΔ)
    plotL(L, T_for_LD, D_for_L, χ, D[d], "Susc", χΔ)
end







#     ---------     Convergence test     ---------     #

E0=load("Data/$(L)_$N/70_0.55_-1.0.jld2")[:"Energies"]
rhat(Chains(E0, ["my param"]))







Gaussmax(T_for_LD[100, 0.01f0],C[7,:,13])
T_for_LD[100,-.035f0][argmax(C[9,:,9])]
T_for_LD[100,.035f0][argmax(χ[9,:,9])]


# @save "AllData.jld2" E EΔ C CΔ M MΔ Mx MxΔ My MyΔ Mz MzΔ χ χΔ χx χxΔ χy χyΔ χz χzΔ Corr T_for_LD D_for_L L