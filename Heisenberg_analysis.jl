using Distributions, Plots, Makie, CairoMakie, ColorSchemes, Colors, CSV, DataFrames, Dierckx, LsqFit, JLD2, Glob, Primes, MCMCChains   # To have basic math (mean), To plot, To plot, To plot (needed for Makie), colors for lattices, write/read csv (CSV & DataFrames), interpolate & fit, save data in compact file, glob to read files and name files, To get divisors of a number, To check convergence
include("./Heisenberg_analysis_functions.jl")


# ----- Initial parameters ----- #


const N    = 1000_000     # Number of lattice sweeps
const L    = [8,12,20,32,50,70]  # Lattices size
const Nperbin = 100       # Number of measurements per bins
const PBC  = true         # Periodic Boundary Conditions
folder_name_end=""    # if the folder as a name like [8, 12, 20]_1000000_Ising, the name end is _Ising. By default the folders are created without name end


# ----- Get data ----- #

# all quantities are intensive (E and M are saved normalized, C and χ are then computed directly normalized)

const T, D, T_for_DL = Get_T_DL(N, L, folder_name_end)   # Obtain Temperature and anisotropic term in Data/$(L)_$N i.e. Data/[8, 12, 20]_1000000
Nmeasurement = length(load("Data/$(L)_$N$(folder_name_end)/$(L[1])_$(T_for_DL[(L[1],D[1])][1])_$(D[1]).jld2")[:"Energies"])
Nbin = div(Nmeasurement, Nperbin)

nL, nT, nD  =   length(L), length(T), length(D) # number of different lattice sizes, temperatures and anisotropic term
E, EΔ       =   zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD)
M, MΔ       =   zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD)
Mx, My, Mz, MxΔ, MyΔ, MzΔ = zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD)
χ, χΔ       =   zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD)
χx, χy, χz, χxΔ, χyΔ, χzΔ       =   zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD)
C, CΔ, CΔ2  =   zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD), zeros(Float32, nL, nT, nD)
Acceptance  =   zeros(4, nL, nT, nD)
Corr        =   Array{Vector{Float32}}(undef, nL, nT, nD)
# AllLattices =   Array{Vector{Array{Float32, 3}}}(undef, nL, nT, nD) # for Correlation time (do not exists by default)
SwapAccept  =   Array{Vector}(undef, nL, nD)    # acceptance of the lattices swap (parallel tempering)

for l in eachindex(L)
    for d in eachindex(D)
        for t in eachindex(T_for_DL[(L[l], D[d])])
            Data = load("Data/$(L)_$N$(folder_name_end)/$(L[l])_$(T_for_DL[(L[l], D[d])][t])_$(D[d]).jld2")
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
            χ[l,t,d], χΔ[l,t,d] = (mean(Ma.^2)-M[l,t,d]^2)/T[t]*L[l]^2, std(Binor2nd(Ma, Nbin, Nperbin, T[t], L[l]))/sqrt(Nbin)
            χx[l,t,d], χxΔ[l,t,d] = (mean(mx.^2)-Mx[l,t,d]^2)/T[t]*L[l]^2, std(Binor2nd(mx, Nbin, Nperbin, T[t], L[l]))/sqrt(Nbin)
            χy[l,t,d], χyΔ[l,t,d] = (mean(my.^2)-My[l,t,d]^2)/T[t]*L[l]^2, std(Binor2nd(my, Nbin, Nperbin, T[t], L[l]))/sqrt(Nbin)
            χz[l,t,d], χzΔ[l,t,d] = (mean(mz.^2)-Mz[l,t,d]^2)/T[t]*L[l]^2, std(Binor2nd(mz, Nbin, Nperbin, T[t], L[l]))/sqrt(Nbin)
            C[l,t,d], CΔ[l,t,d], CΔ2[l,t,d] = (mean(En.^2)-E[l,t,d]^2)/T[t]^2*L[l]^2, Errorpropagation(Binor(En, Nbin, Nperbin), EΔ[l,t,d])/T[t]^2*L[l]^2,      std(Binor2nd(En, Nbin, Nperbin, T[t], L[l]))/T[t]/sqrt(Nbin)
            Corr[l,t,d]         = Data[:"corr"]
            Acceptance[:,l,t,d] = Data[:"accept"]
            # AllLattices[l,t,d] = Data[:"Lattices"]
            # println("N = ", L[l], "\tT = ", T[t], " \td = ", D[d], "\t E = ", round(E[l,t,d];digits=3), " ± ", round(EΔ[l,t,d];digits=5), "\tM = ", round(M[l,t,d];digits=3), " ± ", round(MΔ[l,t,d];digits=5), " \t χ = ", round(χ[l,t,d];digits=3), " ± ", round(χΔ[l,t,d];digits=3), "\tC = ", round(C[l,t,d];digits=3), " ± ", round(CΔ[l,t,d];digits=4), "\taccept = ", round.(Acceptance[:,l,t,d];digits=3))
        end
        SwapAccept[l,d] = load("Data/$(L)_$N$(folder_name_end)/swap_$(L[l])_$(D[d])_$(T[1]).jld2")[:"SwapAccept"]
    end
    print(L[l], "     ")
end


# -----     Plots and analysis     ----- #


PlotColors = cgrad([RGB(.8,.9,1), RGB(0,0,.5), RGB(0,0,0), RGB(.6,0,0), RGB(1,.75,.75)], [0., .4999, .5, .5001, 1.], categorical = false) # Color for plots
colorbar(PlotColors, D) # Need to create the colorbar for some following graphs    

plotL(L, T_for_DL, D, M, 1, "M", MΔ)  # Plot all lattice sizes for the same d.    Fifth entry for which d
plotL(L, T_for_DL, D, Mx, 1, "M_x", MxΔ)
plotL(L, T_for_DL, D, χ, 1, "Susc", χΔ)
plotL(L, T_for_DL, D, χx, 1, "Susc_x", χxΔ)
plotL(L, T_for_DL, D, χy, 1, "Susc_y", χyΔ)
plotL(L, T_for_DL, D, χz, 1, "Susc_z", χzΔ)
plotL(L, T_for_DL, D, E, 1, "E", EΔ)
plotL(L, T_for_DL, D, C, 1, "C", CΔ2)
println("80 swap : $(round.(mean(SwapAccept[:,1]); digits=2))")
plot(T[2:end],diff(E[1,:,1]))

plotD(PlotColors, L, T_for_DL, D, χ, "Susc", χΔ)   # Plot all the d for a given lattice size (last number for lattice size : L[#], the bigger one by default)
# plotD(L, T_for_DL, D, C, "C", CΔ, 3)   # with error propag from ΔE, smaller error bars
plotD(PlotColors, L, T_for_DL, D, C, "C", CΔ2)
plotD(PlotColors, L, T_for_DL, D, E, "E", EΔ, 1)
plotD(PlotColors, L, T_for_DL, D, M, "M", MΔ)


# α, σα, γ, σγ = zeros(nD), zeros(nD), zeros(nD), zeros(nD)
# for i=1:nD
#     α[i], σα[i] = crit(L, T, C[:,:,i])#, "Capacity")
#     γ[i], σγ[i] = crit(L, T, χ[:,:,i])#, "Susceptibility")
# end
# println(α, "\n", σα, "\n\n", γ, "\n", σγ)



ξ, σξ = critlength(T, Corr[1,:,1], 1.0, true)           # third entry for the wanted temperature, fourth to plot or not, fifth component to "true" to get exponent of the algebraic decay instead of correlation length
Plot_Max_ξ(L, D, T_for_DL, Corr, 7, true) # bad for now as it always tries to fit an exponential, but below Tc, this is a power law

# CorrTimePlot(AllLattices, T_for_DL, D, 3, L, 2, false)   # Works only when lattices have been saved in the simulation process, doesn't by default.   Fourth element for which d, sixth for which L

histogram([x for cell in SwapAccept for x in cell])


crit(L, T_for_DL[(L[1], D[1])], C[:,:,1])

interpmax(T, C[2,:,1])

Tmax_χ, χmax, χfitln, χfitpower, χfit_Tc = Plot_Max_C_Χ(PlotColors, D, T_for_DL, χ, χΔ, L, "Susc") # Plot and fit the max of C or χ
Tmax_C, Cmax, Cfitln, Cfitpower, Cfit_Tc = Plot_Max_C_Χ(PlotColors, D, T_for_DL, C, CΔ, L, "C")

a = CritLength(Corr[5,20,13], L[5])

for i in [1, Int((length(D)+1)/2), length(D)]
    println("$(D[i]) χ\tln: $(χfitln[i])\tpower (p2*x^p1), then error: $(χfitpower[i])\t Tmax-Tc propto L^ $(χfit_Tc[i][1][1]) ± $(χfit_Tc[i][2][1])\t Tc = $(χfit_Tc[i][1][2]) ± $(χfit_Tc[i][2][2])")
end
for i in [1, Int((length(D)+1)/2), length(D)]
    println("$(D[i]) C\tln: $(Cfitln[i])\tpower (p2*x^p1), then error: $(Cfitpower[i])")#\t Tmax-Tc propto L^ $(Cfit_Tc[i][1][1]) ± $(Cfit_Tc[i][2][1])\t Tc = $(Cfit_Tc[i][1][2]) ± $(Cfit_Tc[i][2][2])")
end

Plots.plot(D, [x[1][1] for x in χfitpower], seriestype=:scatter, yerr=[x[2][1] for x in χfitpower], label="", title="χ ∝ L^-γ vs D", xlabel="⬅ Ising                              Value of D                              XY ➡", ylabel="-γ")
Plots.plot(D, [x[1][1] for x in Cfitpower], seriestype=:scatter, yerr=[x[2][1] for x in Cfitpower], label="", title="C ∝ L^a vs D", xlabel="⬅ Ising                              Value of D                              XY ➡", ylabel="a")

Plots.plot(D, [x[1] for x in χfitln], seriestype=:scatter, yerr=[x[2] for x in χfitln], label="", title="χ = a*ln(L) vs D", xlabel="⬅ Ising                              Value of D                              XY ➡", ylabel="a")
Plots.plot(D, [x[1] for x in Cfitln], seriestype=:scatter, yerr=[x[2] for x in Cfitln], label="", title="C = a*ln(L) vs D", xlabel="⬅ Ising                              Value of D                              XY ➡", ylabel="a")

Plots.plot(D, [x[1][2] for x in χfit_Tc], seriestype=:scatter, yerr=[x[2][2] for x in χfit_Tc], label="", title = "T of χ_max for L → ∞ vs D", xlabel="⬅ Ising                              Value of D                              XY ➡", ylabel="T of χ_max")


# Plots.savefig("Plot/Susc_Tc.pdf")


#     ---------     Convergence test     ---------     #

E0=load("Data/$(L)_$N/70_0.55_-1.0.jld2")[:"Energies"]
chn = Chains(E0, ["my param"])
rhat(chn)


# fit corrlength (see where critical length diverges)
# correlation vs distance is not correlation length
# Yeomans 116 (115-118) to rescale C and 


for i=1:2:23
    plotL(L, T_for_DL, D, χ, i, "Susc", χΔ)
end




for i=1:23
    println("$i  ",round.(SwapAccept[6,i];digits=3), "\n")
end

a=6
println(round.(Tmax_C[:,a] .- Tmax_χ[:,a];digits=3))
Tmax_C[1,2]

SwapAccept[6,23][10:19]

for d=1:23
    p=Plots.plot(T_for_DL[L[6],D[d]], C[6,1:length(T_for_DL[L[6],D[d]]),d]/maximum(C[6,:,d]))
    Plots.plot!(T_for_DL[L[6],D[d]], χ[6,1:length(T_for_DL[L[6],D[d]]),d]/maximum(χ[6,:,d]))
    display(p)
end

for i=1:23
    println("$(D[i]) \t $(round((Tmax_C[i,6]+2*Tmax_χ[i,6])/3;digits=3))\n")
end