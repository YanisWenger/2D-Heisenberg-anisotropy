using Plots, ColorSchemes, Colors, CSV, DataFrames, Dierckx, LsqFit, JLD2, Glob   # To plot, colors for lattices, write/read csv (CSV & DataFrames), interpolate, save data in compact file, glob to read files and name files
include("./Heisenberg_analysis_functions.jl")

Nbin = 50       # Number of bins
PBC=true        # Periodic Boundary Conditions
N = 1040_000    # Ising
N = 1100_000    # XY
N = 1000_000    # XYZ
L=[8,12,20,32]  # Lattices size
L = [8, 12, 20]

N=20000
L=[8,12]

N = 50_000_000
L=[10]

T, d, T_for_d = Get_T_d(N, L)   # Obtain Temperature and anisotropic term

nL, nT, nd  =   length(L), length(T), length(d) # number of different lattices, temperatures and d
E, EΔ       =   zeros(Float32, nL, nT, nd), zeros(Float32, nL, nT, nd)
M, MΔ       =   zeros(Float32, nL, nT, nd), zeros(Float32, nL, nT, nd)
χ, χΔ       =   zeros(Float32, nL, nT, nd), zeros(Float32, nL, nT, nd)
C, CΔ       =   zeros(Float32, nL, nT, nd), zeros(Float32, nL, nT, nd)
Acceptance  =   zeros(4, nL, nT, nd)
Corr        =   Array{Any}(undef, nL, nT, nd)
# AllLattices =   Array{Any}(undef, nL, nT, nd)
SwapAccept  =   Array{Vector}(undef, nL, nd)
# ----------------------------------------------- #

for l in eachindex(L)
    for z in eachindex(d)
        for t in eachindex(T_for_d[d[z]])
            Data = load("Data/$(L)_$N/$(L[l])_$(T_for_d[d[z]][t])_$(d[z]).jld2")
            E[l,t,z], EΔ[l,t,z] = mean(Data[:"Energies"]),                            std(Binor(Data[:"Energies"], Nbin))
            M[l,t,z], MΔ[l,t,z] = mean(Data[:"Mag"]),                                 std(Binor(Data[:"Mag"], Nbin))
            χ[l,t,z], χΔ[l,t,z] = (mean(Data[:"Mag"].^2)-M[l,t,z]^2)/T[t]*L[l]^2,   std(Binor2nd(Data[:"Mag"], Nbin, T[t], L[l]))
            C[l,t,z], CΔ[l,t,z] = (mean(Data[:"Energies"].^2)-E[l,t,z]^2)/T[t]^2*L[l]^2, Errorpropagation(Binor(Data[:"Energies"], Nbin), EΔ[l,t,z])/T[t]^2*L[l]^2
            Corr[l,t,z]       = Data[:"corr"]
            Acceptance[:,l,t,z] = Data[:"accept"]
            # AllLattices[l,t,z] = Data[:"Lattices"]
            println("N = ", L[l], "\tT = ", T[t], "\td = ", d[z], "\t E = ", round(E[l,t,z];digits=3), " ± ", round(EΔ[l,t,z];digits=3), "\tM = ", round(M[l,t,z];digits=3), " ± ", round(MΔ[l,t,z];digits=3), " \t χ = ", round(χ[l,t,z];digits=7), " ± ", round(χΔ[l,t,z];digits=7), "\tC = ", round(C[l,t,z];digits=2), " ± ", round(CΔ[l,t,z];digits=5), "\taccept = ", round.(Acceptance[:,l,t,z];digits=3))
        end
        SwapAccept[l,z] = load("Data/$(L)_$N/swap_$(L[l])_$(d[z]).jld2")[:"SwapAccept"]
    end
end
println(SwapAccept)

basicplot1(L, T_for_d, d, 1000.0, M, "M", MΔ) # fourth number for the value of d
basicplot1(L, T_for_d, d, 1000.0, χ, "Susc", χΔ)
basicplot1(L, T_for_d, d, 0.0, E, "E", EΔ)
basicplot1(L, T_for_d, d, 0.0, E2, "E2", EΔ2)
basicplot1(L, T_for_d, d, 1000.0, C, "C", CΔ)

basicplotd(T_for_d, d, χ, "Susc", χΔ,1)
basicplotd(T_for_d, d, C, "C", CΔ, 1) # last number for lattice size : L[#]
basicplotd(T_for_d, d, E, "E", EΔ)
basicplotd(T_for_d, d, M, "M", MΔ)


α, σα, γ, σγ = zeros(nd), zeros(nd), zeros(nd), zeros(nd)
for i=1:nd
    α[i], σα[i] = crit(L, T, C[:,:,i])#, "Capacity")
    γ[i], σγ[i] = crit(L, T, χ[:,:,i])#, "Susceptibility")
end
println(α, "\n", σα, "\n\n", γ, "\n", σγ)

ξ, σξ = critlength(T, Corr, 1.0, false) # shouldn't work for now


CorrPlot(AllLattices, T_for_d, d, 3, L, 2, false) # first number for which d, second for which L