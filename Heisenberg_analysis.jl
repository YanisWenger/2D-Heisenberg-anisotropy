using Distributions, Plots, Makie, CairoMakie, ColorSchemes, Colors, CSV, DataFrames, Dierckx, LsqFit, JLD2, Glob, MCMCChains   # To have basic math (mean), To plot, To plot (colorbar for D), To plot (needed for Makie), colors for lattices, write/read csv (CSV & DataFrames), interpolate & fit, save data in compact file, glob to read files and name files, To check convergence
include("./Heisenberg_analysis_functions.jl")


# ----- Initial parameters ----- #


const L = [8,12,20,32,50,70,100,140,200,300]  # Lattices size
const N = 1000_000
const FolderName = "$(L)_$(N)"
const Nperbin = 100       # Number of measurements per bins
const PBC  = true         # Periodic Boundary Conditions

# ----- Get data ----- #

# all quantities are intensive (E and M are saved normalized (that is why we multiply (instead of divide) by L^2 to compute C and χ))

const T, D_for_L, T_for_LD = Get_T_LD(L, FolderName)   # Obtain Temperature and anisotropic term in Data/$(L)_$N
nL, nT, nD  =   length(L), length(T), maximum([length(D_for_L[i]) for i in L]) # number of different lattice sizes, temperatures and anisotropic term
E, EΔ       =   zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD)
M, MΔ       =   zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD)
Mxy, Mz, MxyΔ, MzΔ = zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD)
χ, χΔ       =   zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD)
χxy, χz, χxyΔ, χzΔ  =   zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD)
C, CΔ       =   zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD)
Acceptance  =   zeros(4, nL, nT, nD)
Corr        =   Array{Vector{Float32}}(undef, nL, nT, nD)
CorrFunc    =   Array{Vector{Float32}}(undef, nL, nT, nD)
# AllLattices =   Array{Vector{Array{Float32, 3}}}(undef, nL, nT, nD) # for Correlation time (do not exists by default)
SwapAccept  =   Array{Vector}(undef, nL, nD)    # acceptance of the lattices swap (parallel tempering)

for l in eachindex(L)
    ll = L[l]
    for d in eachindex(D_for_L[ll])
        dd = D_for_L[ll][d]
        for t in eachindex(T_for_LD[(ll, dd)])
            tt = T_for_LD[(ll, dd)][t]
            l2t = ll^2/tt
            Data = load("Data/$FolderName/$(ll)_$(tt)_$dd.jld2")
            En = Data[:"Energies"]
            mz = abs.(Data[:"Mz"])
            mxy = sqrt.(Data[:"Mx"].^2+Data[:"My"].^2)
            Ma = sqrt.(mxy.^2+mz.^2)
            MeanE = mean(En);   MeanMag = mean(Ma)
            MeanMxy = mean(mxy);  MeanMz = mean(mz)
            Nbin = div(length(En), Nperbin)
            E[l,t,d],   EΔ[l,t,d]   = MeanE,                        std(Binor(En, Nbin, Nperbin))/sqrt(Nbin)
            M[l,t,d],   MΔ[l,t,d]   = MeanMag,                      std(Binor(Ma, Nbin, Nperbin))/sqrt(Nbin)
            Mxy[l,t,d], MxyΔ[l,t,d] = MeanMxy,                      std(Binor(mxy, Nbin, Nperbin))/sqrt(Nbin)
            Mz[l,t,d],  MzΔ[l,t,d]  = MeanMz,                       std(Binor(mz, Nbin, Nperbin))/sqrt(Nbin)
            χ[l,t,d],   χΔ[l,t,d]   = (mean(Ma.^2)-MeanMag^2)*l2t,  std(Binor2nd(Ma, Nbin, Nperbin, tt, ll))/sqrt(Nbin) #   multiplied by L^2 as M is normalized
            χxy[l,t,d], χxyΔ[l,t,d] = (mean(mxy.^2)-MeanMxy^2)*l2t, std(Binor2nd(mxy, Nbin, Nperbin, tt, ll))/sqrt(Nbin)
            χz[l,t,d],  χzΔ[l,t,d]  = (mean(mz.^2)-MeanMz^2)*l2t,   std(Binor2nd(mz, Nbin, Nperbin, tt, ll))/sqrt(Nbin)
            C[l,t,d],   CΔ[l,t,d]   = mean((En .-MeanE).^2)*l2t/tt, std(Binor2nd(En.-min(dd, 0), Nbin, Nperbin, tt, ll))/tt/sqrt(Nbin)    #   multiplied by L^2 as C is normalized
            DataCorr                = Data[:"corr"]
            Corr[l,t,d]             = DataCorr
            CorrFunc[l,t,d]         = DataCorr .- MeanMag^2
            Acceptance[:,l,t,d]     = Data[:"accept"]
            # AllLattices[l,t,d] = Data[:"Lattices"]
            # println("N = ", ll, "\tT = ", tt, " \td = ", dd, "\t E = ", round(E[l,t,d];digits=3), " ± ", round(EΔ[l,t,d];digits=5), "\tM = ", round(M[l,t,d];digits=3), " ± ", round(MΔ[l,t,d];digits=5), " \t χ = ", round(χ[l,t,d];digits=3), " ± ", round(χΔ[l,t,d];digits=3), "\tC = ", round(C[l,t,d];digits=3), " ± ", round(CΔ[l,t,d];digits=4), "\taccept = ", round.(Acceptance[:,l,t,d];digits=3))
        end
        # SwapAccept[l,d] = load("Data/$FolderName/swap_$(ll)_$(dd)_$(T_for_LD[(ll, dd)][1]).jld2")[:"SwapAccept"]
    end
    print(ll, "     ")
end


# X = load("AllDataAugust.jld2")
# E=X["E"]; EΔ=X["EΔ"]; C=X["C"]; CΔ=X["CΔ"]; M=X["M"]; MΔ=X["MΔ"]; Mxy=X["Mxy"]; MxyΔ=X["MxyΔ"]; Mz=X["Mz"]; MzΔ=X["MzΔ"]; χ=X["χ"]; χΔ=X["χΔ"]; χxy=X["χxy"]; χxyΔ=X["χxyΔ"]; χz=X["χz"]; χzΔ=X["χzΔ"]; Corr=X["Corr"]; const T_for_LD=X["T_for_LD"]; const L=X["L"]; const D_for_L=X["D_for_L"]; CorrFunc=X["CorrFunc"] 



# -----     Plots and analysis     ----- #



default(titlefont = (18, "Computer Modern"), xguidefont = (14, "Computer Modern"), yguidefont = (14, "Computer Modern"), xtickfont = (12, "Computer Modern"), ytickfont = (12, "Computer Modern"), legendfont=(12, "Computer Modern"))
PlotColors = cgrad([RGB(.5,.9,1), RGB(0,0,.5), RGB(0,0,0), RGB(.6,0,0), RGB(1,.75,.75)], [0., .4999, .5, .5001, 1.], categorical = false) # Color for plots
colorbar(PlotColors, D_for_L) # Need to create the colorbar for some following graphs    

plotL(L, T_for_LD, D_for_L, M, 0, "Magnetization", MΔ)  # Plot all lattice sizes for the same d.    Fifth entry for which d
plotL(L, T_for_LD, D_for_L, Mz, -1000, "Magnetization along \$\\hat{z}\$", MzΔ)
plotL(L, T_for_LD, D_for_L, Mxy, 1000, "Magnetization on XY plane", MxyΔ)
plotL(L, T_for_LD, D_for_L, χ, 0f0, "Susceptibility", χΔ)
plotL(L, T_for_LD, D_for_L, χxy, 0, "Susceptibility on XY plane", χxyΔ)
plotL(L, T_for_LD, D_for_L, χz, 0, "Susceptibility along \$\\hat{z}\$", χzΔ)
plotL(L, T_for_LD, D_for_L, E, 0, "Energy", EΔ)
plotL(L, T_for_LD, D_for_L, C, 0f0, "Heat capacity", CΔ)


plotD(PlotColors, L, T_for_LD, D_for_L, χ, "Susc", χΔ)   # Plot all the d for a given lattice size (last number for lattice size : L[#], the bigger one by default)
# plotD(L, T_for_LD, D, C, "C", CΔ2, 3)   # with error propag from ΔE, smaller error bars
plotD(PlotColors, L, T_for_LD, D_for_L, C, "C", CΔ)
plotD(PlotColors, L, T_for_LD, D_for_L, E, "E", EΔ, 1)
plotD(PlotColors, L, T_for_LD, D_for_L, M, "M", MΔ)

Plots.savefig("Plot/E_Heisenberg.pdf")

plotTCorrelation(L, T_for_LD, D_for_L, Corr, 1000, 4)  # second to last: D value, last: lattice number
plotLCorrelation(L, T_for_LD, D_for_L, Corr, 1000, .8)  # second to last: D value, last: temperature
plotDCorrelation(PlotColors, L, T_for_LD, D_for_L, Corr, 4, .8)  # second to last: lattice number, last: temperature


# crit(L, T_for_LD[(L[1], D[1])], C[:,:,1])
# α, σα, γ, σγ = zeros(nD), zeros(nD), zeros(nD), zeros(nD)
# for i=1:nD
#     α[i], σα[i] = crit(L, T, C[:,:,i])#, "Capacity")
#     γ[i], σγ[i] = crit(L, T, χ[:,:,i])#, "Susceptibility")
# end
# println(α, "\n", σα, "\n\n", γ, "\n", σγ)



ξ, σξ, p = critlength(Corr, T_for_LD, L, D_for_L, 9, 0.7, 1.0, true, false)           # sixth entry for the wanted temperature, seventh to plot or not, eighth component to "true" to get exponent of the algebraic decay instead of correlation length

let
l=7; d=-1
ξ, σξ, p1 = critlength(CorrFunc, T_for_LD, L, D_for_L, l, d, 0.0, false, false)           # sixth entry for the wanted temperature, seventh to plot or not, eighth component to "true" to get exponent of the algebraic decay instead of correlation length
ξ, σξ, p2 = critlength(CorrFunc, T_for_LD, L, D_for_L, l, d, 0.0, false, true)           # sixth entry for the wanted temperature, seventh to plot or not, eighth component to "true" to get exponent of the algebraic decay instead of correlation length
ξ, σξ, p3 = critlength(CorrFunc, T_for_LD, L, D_for_L, l, d, 10, false, false)           # sixth entry for the wanted temperature, seventh to plot or not, eighth component to "true" to get exponent of the algebraic decay instead of correlation length
ξ, σξ, p4 = critlength(CorrFunc, T_for_LD, L, D_for_L, l, d, 10, false, true)           # sixth entry for the wanted temperature, seventh to plot or not, eighth component to "true" to get exponent of the algebraic decay instead of correlation length
p=Plots.plot(p1, p2, p3, p4, layout=(2, 2), size=(600,400))#, size=(800, 600), sgtitle="Correlation for L=$(L[l]), D=$(d)")
display(p)
# Plots.savefig("Plot/Correlation_function_Ising.pdf")
# Plots.savefig("Plot/Correlation_$(L[l])_$(d).pdf")
end

Plots.plot(CorrFunc[7,20,1], label="T=$(T_for_LD[100, 1f3][20])", xlabel = "Distance", ylabel="Correlation function", framestyle = :box)
Plots.savefig("Plot/Correlation_function_XY_100_1.008.pdf")

# Plot_Max_ξ(L, D, T_for_LD, Corr, 7, true) # bad for now as it always tries to fit an exponential, but below Tc, this is a power law

# CorrTimePlot(AllLattices, T_for_LD, D, 3, L, 2, false)   # Works only when lattices have been saved in the simulation process, doesn't by default.   Fourth element for which d, sixth for which L

histogram([x for cell in SwapAccept for x in cell])


Tmax_χ, χmax, χfitln, χfitln_err, χfitpower, χfit_Tc = Plot_Max_C_Χ(PlotColors, D_for_L, T_for_LD, χ, χΔ, L, "Susc") # Plot and fit the max of C or χ
Tmax_C, Cmax, Cfitln, Cfitln_err, Cfitpower, Cfit_Tc = Plot_Max_C_Χ(PlotColors, D_for_L, T_for_LD, C, CΔ, L, "C")


Plots.plot(D, [x[1][1] for x in χfitpower], seriestype=:scatter, yerr=[x[2][1] for x in χfitpower], label="", title="χ ∝ L^-γ vs D", framestyle = :box, xlabel="⬅ Ising                              Value of D                              XY ➡", ylabel="-γ")
Plots.plot(D, [x[1][1] for x in Cfitpower], seriestype=:scatter, yerr=[x[2][1] for x in Cfitpower], label="", title="C ∝ L^a vs D", framestyle = :box, xlabel="⬅ Ising                              Value of D                              XY ➡", ylabel="a")

Plots.plot(D, [x[1] for x in χfitln], seriestype=:scatter, yerr=[x[2] for x in χfitln], label="", title="χ = a*ln(L) vs D", framestyle = :box, xlabel="⬅ Ising                              Value of D                              XY ➡", ylabel="a")
Plots.plot(D, Cfitln[:,1], seriestype=:scatter, yerr=Cfitln_err[:,1], label="", title="C = a*ln(L/L0)*(1+L'/L) vs D", framestyle = :box, xlabel="⬅ Ising                              Value of D                              XY ➡", ylabel="a")
Plots.plot(D, Cfitln[:,2], seriestype=:scatter, yerr=Cfitln_err[:,2], label="", title="C = a*ln(L/L0)*(1+L'/L) vs D", framestyle = :box, xlabel="⬅ Ising                              Value of D                              XY ➡", ylabel="L0")
Plots.plot(D, Cfitln[:,3], seriestype=:scatter, yerr=Cfitln_err[:,3], label="", title="C = a*ln(L/L0)*(1+L'/L) vs D", framestyle = :box, xlabel="⬅ Ising                              Value of D                              XY ➡", ylabel="L'")

Plots.plot(D, [x[1][2] for x in χfit_Tc], seriestype=:scatter, yerr=[x[2][2] for x in χfit_Tc], label="", title = "T of χ_max for L → ∞ vs D", framestyle = :box, xlabel="⬅ Ising                              Value of D                              XY ➡", ylabel="T of χ_max")


scaling_plot_C_Ising(L, T_for_LD, D_for_L, C, -0.7, Tmax_C, CΔ)
scaling_plot_χ_Ising(L, T_for_LD, D_for_L, χ, -0.7, Tmax_χ, χΔ)
scaling_plot_C_XY(L, T_for_LD, D_for_L, C, 1000, Tmax_C, CΔ)
scaling_plot_χ_XY(L, T_for_LD, D_for_L, χ, 1000, Tmax_χ, χΔ)
Plots.savefig("Plot/C_aln(L).pdf")


# fit corrlength (see where critical length diverges)




# @save "AllDataAugust.jld2" E EΔ C CΔ M MΔ Mxy MxyΔ Mz MzΔ χ χΔ χxy χxyΔ χz χzΔ Corr T_for_LD D_for_L L CorrFunc



LogLogPlotAndFit(L, [χmax[l,0] for l in L], [1.75, .7], "Lattice length", "Maximum of Susceptibility","\$\\chi_{max}\$")
LogLogPlotAndFit(L, [Cmax[l,0] for l in L], [1.75, .7], "Lattice length", "Maximum of Heat capacity","\$C_{max}\$")
xlogPlotAndFit(L, [Cmax[l,-1000] for l in L], [1.75, .7])

Plots.savefig("Plot/Susc_max_d0.pdf")







#     ---------     Convergence test     ---------     #





rhat(Chains(load("Data/[400]_1000000_new/400_0.631_-0.02.jld2")["Energies"], ["my param"])).nt[2][1] # Gelman-Rubin convergencee test


conv = Convergence_check(L, "$(L)_$(N)", 300)  # all lattice sizes in the folder, folder name, lattice size that we want to know the convergence (optional, if none, calculate for all lattice sizes)