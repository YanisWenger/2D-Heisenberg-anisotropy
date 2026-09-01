using Distributions, Plots, Makie, CairoMakie, ColorSchemes, Colors, CSV, DataFrames, Dierckx, LsqFit, JLD2, Glob, MCMCChains   # To have basic math (mean), To plot, To plot (colorbar for D), To plot (needed for Makie), colors for lattices, write/read csv (CSV & DataFrames), interpolate & fit, save data in compact file, glob to read files and name files, To check convergence
include("./Heisenberg_analysis_functions.jl")


# ----- Initial parameters ----- #


const L = [8,12,20,32,50,70,100]#,140,200,300,400]  # Lattices size
const N = 1000_000
const FolderName = "$(L)_$(N)_Ising"
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
ConnCorr    =   Array{Vector{Float32}}(undef, nL, nT, nD)
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
            ConnCorr[l,t,d]         = DataCorr .- MeanMag^2
            Acceptance[:,l,t,d]     = Data[:"accept"]
            # AllLattices[l,t,d] = Data[:"Lattices"]
        end
        # SwapAccept[l,d] = load("Data/$FolderName/swap_$(ll)_$(dd)_$(T_for_LD[(ll, dd)][1]).jld2")[:"SwapAccept"]  # acceptance for swapping configuration in the parallel tempering algorithm
    end
    print(ll, "     ")
end

# @save "AllData.jld2" E EΔ C CΔ M MΔ Mxy MxyΔ Mz MzΔ χ χΔ χxy χxyΔ χz χzΔ Corr T_for_LD D_for_L L ConnCorr # save the computed data to save some time next time


# X = load("AllDataAugust.jld2")    # get the data that was already saved
# E=X["E"]; EΔ=X["EΔ"]; C=X["C"]; CΔ=X["CΔ"]; M=X["M"]; MΔ=X["MΔ"]; Mxy=X["Mxy"]; MxyΔ=X["MxyΔ"]; Mz=X["Mz"]; MzΔ=X["MzΔ"]; χ=X["χ"]; χΔ=X["χΔ"]; χxy=X["χxy"]; χxyΔ=X["χxyΔ"]; χz=X["χz"]; χzΔ=X["χzΔ"]; Corr=X["Corr"]; const T_for_LD=X["T_for_LD"]; const L=X["L"]; const D_for_L=X["D_for_L"]; ConnCorr=X["ConnCorr"] 



# -----     Plots and analysis     ----- #



default(titlefont = (18, "Computer Modern"), xguidefont = (14, "Computer Modern"), yguidefont = (14, "Computer Modern"), xtickfont = (12, "Computer Modern"), ytickfont = (12, "Computer Modern"), legendfont=(12, "Computer Modern"), framestyle = :box)
PlotColors = cgrad([RGB(.5,.9,1), RGB(0,0,.5), RGB(0,0,0), RGB(.6,0,0), RGB(1,.75,.75)], [0., .4999, .5, .5001, 1.], categorical = false) # Color for plots
colorbar(PlotColors, D_for_L) # Need to create the colorbar for some following graphs    

plotL(L, T_for_LD, D_for_L, M, 0, "Magnetization", MΔ)  # Plot all lattice sizes for the same d.    Fifth entry for which d
plotL(L, T_for_LD, D_for_L, Mz, -1000, "Magnetization along \$\\hat{z}\$", MzΔ)
plotL(L, T_for_LD, D_for_L, Mxy, 1000, "Magnetization on XY plane", MxyΔ)
plotL(L, T_for_LD, D_for_L, χ, -.02f0, "Susceptibility", χΔ)
plotL(L, T_for_LD, D_for_L, χxy, -.02f0, "Susceptibility on XY plane", χxyΔ)
plotL(L, T_for_LD, D_for_L, χz, -.02f0, "Susceptibility along \$\\hat{z}\$", χzΔ)
plotL(L, T_for_LD, D_for_L, E, 0, "Energy", EΔ)
plotL(L, T_for_LD, D_for_L, C, 0f0, "Heat capacity", CΔ)


plotD(PlotColors, L, T_for_LD, D_for_L, χ, "Susceptibility", χΔ, 9)   # Plot all the d for a given lattice size (last number for lattice size : L[#], the bigger one by default)
plotD(PlotColors, L, T_for_LD, D_for_L, C, "Heat capacity", CΔ, 9)
plotD(PlotColors, L, T_for_LD, D_for_L, E, "Energy", EΔ, 7)
plotD(PlotColors, L, T_for_LD, D_for_L, M, "Magnetization", MΔ,7)

Plots.savefig("Plot/M_100.pdf")

plotTCorrelation(L, T_for_LD, D_for_L, Corr, .1, 70)  # second to last: D value, last: lattice number
plotLCorrelation(L, T_for_LD, D_for_L, Corr, .1, .8)  # second to last: D value, last: temperature
plotDCorrelation(PlotColors, L, T_for_LD, D_for_L, Corr, 32, .8)  # second to last: lattice number, last: temperature

# d = 23; l=1; ttt = T_for_LD[L[l],D_for_L[L[l]][d]]; p=Plots.plot(ttt, χ[l,1:length(ttt),d], color=:green, label="", xaxis="Temperature", yaxis = "X value")




ξ, σξ = critlength(ConnCorr, T_for_LD, L, D_for_L, 100, -1000f0, 0, false, true, "Connected Correlation function")           # sixth entry for the wanted temperature, seventh to plot or not, eighth component to "true" to get exponent of the algebraic decay instead of correlation length

Plots.plot(ConnCorr[7,20,1], label="T=$(T_for_LD[100, 1f3][20])", xlabel = "Distance", ylabel="Correlation function", framestyle = :box)
Plots.savefig("Plot/Conn_Corr_Ising_low.pdf")


Fit_ξ(L, D_for_L, T_for_LD, Corr, 100, -1000f0)
critic_exp, critic_pow, Tmax_ξ = Fit_ξ(L, D_for_L, T_for_LD, ConnCorr)
# It seems there is always a max for ConnCorr when fitted exp
# "Divergence" in Corr length for Corr (not ConnCorr)

for dd in D_for_L[8]
    aaa = []
    println(dd)
    for ll in L[5:9]
        push!(aaa, Tmax_ξ[ll,dd])
    end
    display(Plots.scatter(L[5:9], aaa, label=dd))#, xaxis=:log, yaxis=:log))
end
# CorrTimePlot(AllLattices, T_for_LD, D, 3, L, 2, false)   # Works only when lattices have been saved in the simulation process, doesn't by default.   Fourth element for which d, sixth for which L

histogram([x for cell in SwapAccept for x in cell])


Tmax_χ, χmax, χfitpower, χfit_Tc = Plot_Max_C_Χ(PlotColors, D_for_L, T_for_LD, χ, χΔ, L, "Susceptibility") # Plot and fit the max of C or χ
Tmax_C, Cmax, Cfitln,    Cfit_Tc = Plot_Max_C_Χ(PlotColors, D_for_L, T_for_LD, C, CΔ, L, "Heat capacity")

if length(values(Cfitln)) == length(D_for_L[L[1]])
    χfitpower_vect=zeros(length(values(χfitpower))); χfitpower_err=zeros(length(values(χfitpower)))
    Cfitln_vect = zeros(length(values(Cfitln)));     Cfitln_err = zeros(length(values(Cfitln)))
    for i = 1:length(D_for_L[L[1]])
        χfitpower_vect[i] = χfitpower[sort(D_for_L[L[1]])[i]][1][1]; χfitpower_err[i] = χfitpower[sort(D_for_L[L[1]])[i]][2][1]
        Cfitln_vect[i] = Cfitln[sort(D_for_L[L[1]])[i]][1][1]; Cfitln_err[i] = Cfitln[sort(D_for_L[L[1]])[i]][2][1]
    end
else
    println("Length of D_for_L[L[1]] ( $(length(D_for_L[L[1]])) ) does not match length of χfitpower ( $(length(values(Cfitln))) )")
end

Plots.plot(D_for_L[L[1]], χfitpower_vect, seriestype=:scatter, yerr=χfitpower_err, label="", title="\$\\chi \\propto L^a\\; vs \\;D\$", framestyle = :box, xlabel="⬅ Ising                Value of D                XY ➡", ylabel="\$a\$")
Plots.plot(D_for_L[L[1]], Cfitln_vect, seriestype=:scatter, yerr=Cfitln_err, label="", title="\$C = a\\cdot \\ln(L/L_0) + L'/L\\;  vs\\; D\$", framestyle = :box, xlabel="⬅ Ising                Value of D                XY ➡", ylabel="\$a\$")


C_Tc_keys = sort([x for x in keys(Cfit_Tc)])
C_Tc_vals = [];     for i in C_Tc_keys;     push!(C_Tc_vals, Cfit_Tc[i]);   end
χ_Tc_keys = sort([x for x in keys(χfit_Tc)])
χ_Tc_vals = [];     for i in χ_Tc_keys;     push!(χ_Tc_vals, χfit_Tc[i]);   end
Plots.plot(χ_Tc_keys, [x[1][2] for x in χ_Tc_vals], seriestype=:scatter, yerr=[x[2][2] for x in χ_Tc_vals], label="\$\\chi\$", title = "Extrapolation of T peak as L → ∞ vs D", framestyle = :box, xlabel="⬅ Ising                Value of D                XY ➡", ylabel="Temperature of the peak")
Plots.plot!(C_Tc_keys, [x[1][2] for x in C_Tc_vals], seriestype=:scatter, yerr=[x[2][2] for x in C_Tc_vals], label="\$C\$")

Plots.plot(χ_Tc_keys, [x[1][1] for x in χ_Tc_vals], seriestype=:scatter, yerr=[x[2][1] for x in χ_Tc_vals], label="\$\\chi\$", title = "Critical exponent", framestyle = :box, xlabel="⬅ Ising                Value of D                XY ➡", ylabel="Critical exponent", ylim = (-2.1,-.1))
Plots.plot!(C_Tc_keys, [x[1][1] for x in C_Tc_vals], seriestype=:scatter, yerr=[x[2][1] for x in C_Tc_vals], label="\$C\$")



scaling_plot_C_Ising(L, T_for_LD, D_for_L, C, -0.7, Tmax_C, CΔ, Cfitln)
scaling_plot_χ_Ising(L, T_for_LD, D_for_L, χ, -0.7, Tmax_χ, χΔ, χfitpower)
scaling_plot_C_XY(   L, T_for_LD, D_for_L, C, 1, Tmax_C, CΔ)
scaling_plot_χ_XY(   L, T_for_LD, D_for_L, χ, 1000, Tmax_χ, χΔ, χfitpower)
Plots.savefig("Plot/Susc_XY_rescaled.pdf")


# fit corrlength (see where critical length diverges)







LogLogPlotAndFit(L, [χmax[l,0] for l in L], [1.75, .7], "Lattice length", "Maximum of susceptibility","\$\\chi_{max}\$")
LogLogPlotAndFit(L, [Cmax[l,0] for l in L], [1.75, .7], "Lattice length", "Maximum of heat capacity","\$C_{max}\$")
xlogPlotAndFit(L, [Cmax[l,-1000] for l in L], [1.75, .7])

Plots.savefig("Plot/Susc_max_d0.pdf")







#     ---------     Convergence test     ---------     #





rhat(Chains(load("Data/[400]_1000000/400_0.631_-0.02.jld2")["Energies"], ["my param"])).nt[2][1] # Gelman-Rubin convergencee test


conv = Convergence_check(L, "$(L)_$(N)", 400)  # 1) all lattice sizes in the folder 2) folder name 3) lattice size that we want to know the convergence (optional, if none, compute for all lattice sizes).     Return in a Vector of tuple the Gelamn-Rubin test value with the corresponding lattice size, temperature and anysotropic term
maximum([x[1] for x in conv])





###     the code assume all value of D are used for the smallest L (L[1])