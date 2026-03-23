using Distributions, Plots, ColorSchemes, Colors, CSV, DataFrames, Dierckx, LsqFit, JLD2, Glob, Primes   # To have basic math (mean), To plot, colors for lattices, write/read csv (CSV & DataFrames), interpolate, save data in compact file, glob to read files and name files, To get divisors of a number
include("./Heisenberg_analysis_functions.jl")

Nbin = 9000     # Number of bins (has to be a divisor of Nmeasurement = (Nstep-burn)/Skip)
PBC=true        # Periodic Boundary Conditions
N = 1040_000    # Ising
N = 1100_000    # XY
N = 1000_000    # XYZ
L = [8,12,20,32]  # Lattices size
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
C, CΔ, CΔ2  =   zeros(Float32, nL, nT, nd), zeros(Float32, nL, nT, nd), zeros(Float32, nL, nT, nd)
Acceptance  =   zeros(4, nL, nT, nd)
Corr        =   Array{Vector{Float32}}(undef, nL, nT, nd)
# AllLattices =   Array{Any}(undef, nL, nT, nd) # for Correlation time (do not exists by default)
SwapAccept  =   Array{Vector}(undef, nL, nd)    # acceptance of the lattices swap (parallel tempering)

# -----  -----   -----    -----     -----    -----   -----  ----- #

for l in eachindex(L)
    for z in eachindex(d)
        for t in eachindex(T_for_d[d[z]])
            Data = load("Data/$(L)_$N/$(L[l])_$(T_for_d[d[z]][t])_$(d[z]).jld2")
            En = Data[:"Energies"]
            Ma = Data[:"Mag"]
            E[l,t,z], EΔ[l,t,z] = mean(En),                             std(Binor(En, Nbin))
            M[l,t,z], MΔ[l,t,z] = mean(Ma),                             std(Binor(Ma, Nbin))
            χ[l,t,z], χΔ[l,t,z] = (mean(Ma.^2)-M[l,t,z]^2)/T[t]*L[l]^2, std(Binor2nd(Ma, Nbin, T[t], L[l]))
            C[l,t,z], CΔ[l,t,z], CΔ2[l,t,z] = (mean(En.^2)-E[l,t,z]^2)/T[t]^2*L[l]^2, Errorpropagation(Binor(En, Nbin), EΔ[l,t,z])/T[t]^2*L[l]^2,      std(Binor2nd(En, Nbin, T[t], L[l]))/T[t]^2
            Corr[l,t,z]         = Data[:"corr"]
            Acceptance[:,l,t,z] = Data[:"accept"]
            # AllLattices[l,t,z] = Data[:"Lattices"]
            println("N = ", L[l], "\tT = ", T[t], " \td = ", d[z], "\t E = ", round(E[l,t,z];digits=3), " ± ", round(EΔ[l,t,z];digits=3), "\tM = ", round(M[l,t,z];digits=3), " ± ", round(MΔ[l,t,z];digits=3), " \t χ = ", round(χ[l,t,z];digits=7), " ± ", round(χΔ[l,t,z];digits=7), "\tC = ", round(C[l,t,z];digits=2), " ± ", round(CΔ[l,t,z];digits=5), "\taccept = ", round.(Acceptance[:,l,t,z];digits=3))
        end
        SwapAccept[l,z] = load("Data/$(L)_$N/swap_$(L[l])_$(d[z]).jld2")[:"SwapAccept"]
    end
end
println(SwapAccept)

plotL(L, T_for_d, d, 1000.0, M, "M", MΔ)    # fourth number for the value of d
plotL(L, T_for_d, d, 1000.0, χ, "Susc", χΔ)
plotL(L, T_for_d, d, 0.0, E, "E", EΔ)
plotL(L, T_for_d, d, 1000.0, C, "C", CΔ)

plotd(L, T_for_d, d, χ, "Susc", χΔ,2)
plotd(L, T_for_d, d, C, "C", CΔ, 3)     # last number for lattice size : L[#]
plotd(L, T_for_d, d, C, "C", CΔ2, 2)    # last number for lattice size : L[#]
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


CorrTimePlot(AllLattices, T_for_d, d, 3, L, 2, false)   # first number for which d, second for which L

swap = [x for cell in SwapAccept for x in cell]
histogram(swap)



crit(L, T_for_d[d[1]], C[:,:,1])

interpmax(2, T, C[:,:,1])

Plot_Max_C_Χ(d, T_for_d, χ, L, "Susc", true)
Plot_Max_C_Χ(d, T_for_d, C, L, "C", true)

Plot_Max_ξ(d, T_for_d, Corr[end,:,:])















function whichbin(vect)
    n = sort(divisors(length(vect)))
    bin = zeros(length(n)-1)
    for i=1:length(bin)
        bin[75-i] = std(Binor2nd(vect, n[i+1],.772f0, 12))
    end
    p=plot(n[1:74],bin,  xaxis = :log10)
    title!("Error bar vs measurement's number per bin)")
    # savefig("Plot/Errorbars.pdf")
    display(p)
    return maximum(bin), argmax(bin)
end

Bestbin = zeros(Float32, nL, 25, nd)
for l in eachindex(L)
    for z in eachindex(d)
        for t=1:25
            Data = load("Data/$(L)_$N/$(L[l])_$(T_for_d[d[z]][t])_$(d[z]).jld2")[:"Energies"]
            a, Bestbin[l,t,z] = whichbin(Data)
        end
    end
end
println(Bestbin)
minimum(Bestbin)
maximum(Bestbin)
mean(Bestbin)
n = sort(divisors(90000))
n[5]

Data = load("Data/$(L)_$N/12_0.65_0.2.jld2")[:"Energies"]
a, Bestbin = whichbin(Data)