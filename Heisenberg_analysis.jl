using Distributions, Plots, ColorSchemes, Colors, CSV, DataFrames, Dierckx, LsqFit, JLD2, Glob, Primes   # To have basic math (mean), To plot, colors for lattices, write/read csv (CSV & DataFrames), interpolate, save data in compact file, glob to read files and name files, To get divisors of a number
include("./Heisenberg_analysis_functions.jl")


# ----- Initial parameters ----- #


Nperbin = 100       # Number of measurements per bins
PBC  = true         # Periodic Boundary Conditions
N    = 1000_000     # Number of lattice sweeps
L    = [8, 12, 20, 32, 50]  # Lattices size


# ----- Get data ----- #

folder_name_end=""
T, d, T_for_d = Get_T_d(N, L, folder_name_end)   # Obtain Temperature and anisotropic term in Data/$(L)_$N i.e. Data/[8, 12, 20]_1000000
Nmeasurement = length(load("Data/$(L)_$N$(folder_name_end)/$(L[1])_$(T_for_d[d[1]][1])_$(d[1]).jld2")[:"Energies"])
Nbin = div(Nmeasurement, Nperbin)

nL, nT, nd  =   length(L), length(T), length(d) # number of different lattices, temperatures and anisotropic term
E, EΔ       =   zeros(Float32, nL, nT, nd), zeros(Float32, nL, nT, nd)
M, MΔ       =   zeros(Float32, nL, nT, nd), zeros(Float32, nL, nT, nd)
χ, χΔ       =   zeros(Float32, nL, nT, nd), zeros(Float32, nL, nT, nd)
C, CΔ, CΔ2  =   zeros(Float32, nL, nT, nd), zeros(Float32, nL, nT, nd), zeros(Float32, nL, nT, nd)
Acceptance  =   zeros(4, nL, nT, nd)
Corr        =   Array{Vector{Float32}}(undef, nL, nT, nd)
# AllLattices =   Array{Vector{Array{Float32, 3}}}(undef, nL, nT, nd) # for Correlation time (do not exists by default)
SwapAccept  =   Array{Vector}(undef, nL, nd)    # acceptance of the lattices swap (parallel tempering)

for l in eachindex(L)
    for z in eachindex(d)
        for t in eachindex(T_for_d[d[z]])
            Data = load("Data/$(L)_$N$(folder_name_end)/$(L[l])_$(T_for_d[d[z]][t])_$(d[z]).jld2")
            En = Data[:"Energies"]
            Ma = Data[:"Mag"]
            E[l,t,z], EΔ[l,t,z] = mean(En),                             std(Binor(En, Nbin, Nperbin))/sqrt(Nbin)
            M[l,t,z], MΔ[l,t,z] = mean(Ma),                             std(Binor(Ma, Nbin, Nperbin))/sqrt(Nbin)
            χ[l,t,z], χΔ[l,t,z] = (mean(Ma.^2)-M[l,t,z]^2)/T[t]*L[l]^2, std(Binor2nd(Ma, Nbin, Nperbin, T[t], L[l]))/sqrt(Nbin)
            C[l,t,z], CΔ[l,t,z], CΔ2[l,t,z] = (mean(En.^2)-E[l,t,z]^2)/T[t]^2*L[l]^2, Errorpropagation(Binor(En, Nbin, Nperbin), EΔ[l,t,z])/T[t]^2*L[l]^2,      std(Binor2nd(En, Nbin, Nperbin, T[t], L[l]))/T[t]/sqrt(Nbin)
            Corr[l,t,z]         = Data[:"corr"]
            Acceptance[:,l,t,z] = Data[:"accept"]
            # AllLattices[l,t,z] = Data[:"Lattices"]
            # println("N = ", L[l], "\tT = ", T[t], " \td = ", d[z], "\t E = ", round(E[l,t,z];digits=3), " ± ", round(EΔ[l,t,z];digits=5), "\tM = ", round(M[l,t,z];digits=3), " ± ", round(MΔ[l,t,z];digits=5), " \t χ = ", round(χ[l,t,z];digits=3), " ± ", round(χΔ[l,t,z];digits=3), "\tC = ", round(C[l,t,z];digits=3), " ± ", round(CΔ[l,t,z];digits=4), "\taccept = ", round.(Acceptance[:,l,t,z];digits=3))
        end
        SwapAccept[l,z] = load("Data/$(L)_$N$(folder_name_end)/swap_$(L[l])_$(d[z]).jld2")[:"SwapAccept"]
    end
    print(L[l], "     ")
end


# -----     Plots and analysis     ----- #



plotL(L, T_for_d, d, M, 1, "M", MΔ)  # Plot all lattice sizes for the same d.    Fifth entry for which d
plotL(L, T_for_d, d, χ, 1, "Susc", χΔ)
plotL(L, T_for_d, d, E, 1, "E", EΔ)
plotL(L, T_for_d, d, C, 1, "C", CΔ2)

plotd(L, T_for_d, d, χ, "Susc", χΔ,5)   # Plot all the d for a given lattice size (last number for lattice size : L[#], the bigger one by default)
# plotd(L, T_for_d, d, C, "C", CΔ, 3)   # with error propag from ΔE, smaller error bars
plotd(L, T_for_d, d, C, "C", CΔ2, 5)
plotd(L, T_for_d, d, E, "E", EΔ)
plotd(L, T_for_d, d, M, "M", MΔ)


α, σα, γ, σγ = zeros(nd), zeros(nd), zeros(nd), zeros(nd)
for i=1:nd
    α[i], σα[i] = crit(L, T, C[:,:,i])#, "Capacity")
    γ[i], σγ[i] = crit(L, T, χ[:,:,i])#, "Susceptibility")
end
println(α, "\n", σα, "\n\n", γ, "\n", σγ)




ξ, σξ = critlength(T, Corr[1,:,1], 1.0, true)           # third entry for the wanted temperature, fourth to plot or not, fifth component to "true" to get exponent of the algebraic decay instead of correlation length
Plot_Max_ξ(d, T_for_d, Corr[1,:,:])

# CorrTimePlot(AllLattices, T_for_d, d, 3, L, 2, false)   # Works only when lattices have been saved in the simulation process, doesn't by default.   Fourth element for which d, sixth for which L

histogram([x for cell in SwapAccept for x in cell])



crit(L, T_for_d[d[1]], C[:,:,1])

interpmax(2, T, C[:,:,1])

Plot_Max_C_Χ(d, T_for_d, χ, L, "Susc")
Plot_Max_C_Χ(d, T_for_d, C, L, "C")

Plot_Max_ξ(d, T_for_d, Corr[end,:,:], false) # bad for now as it always tries to fit an exponential, but below Tc, this is a power law

a = CritLength(Corr[5,20,13], L[5])



# error bar on Plot_Max_C_Χ
# fit corrlength