mutable struct Replica
    Lattice::Array{Float32,3}
    E::Float32          # Instantaneous
    T::Float32
    σ::Float32          # Instantaneous
    Acceptance::Vector{Int64}   # History
    Try::Vector{Int64}          # History
end

function Energy(Lattice::Array{Float32, 3}, L::Int, PBC::Bool, d::Float32)    # Calculate the Energy of a configuration
    E=0
    n = PBC ? L : L-1
    for i=1:n
        for j=1:n
            E += -Correlation(Lattice[:,i,j],Lattice[:,i,mod(j,L)+1]) - Correlation(Lattice[:,i,j],Lattice[:,mod(i,L)+1,j]) + d*cos(Lattice[2,i,j])^2
        end
    end
    return E
end

function Initial_lattice(L::Int, pi32::Float32)
    lattice = Array{Float32}(undef, 2, L, L)
    lattice[1, :, :] .= 2*pi32*rand(TaskLocalRNG(), Float32, L, L)
    lattice[2, :, :] .= acos.(1 .- 2*rand(TaskLocalRNG(), Float32, L, L))
    return lattice
end

function Correlation(spin1::Vector, spin2::Vector)
    return sin(spin1[2])*sin(spin2[2])*cos(spin1[1]-spin2[1]) + cos(spin1[2])*cos(spin2[2])  
end

function EnergyDiff(spin0old::Vector, spin0new::Vector, spin1::Vector, spin2::Vector, spin3::Vector, spin4::Vector, d::Float32) # calculate once sin and cos of old and new spin instead of 4, and once spin and cos of spin2-3-4-5 instead if 2
    sold, snew = sin(spin0old[2]), sin(spin0new[2])
    cold, cnew = cos(spin0old[2]), cos(spin0new[2])
    return (snew*cos(spin0new[1]-spin1[1]) - sold*cos(spin0old[1]-spin1[1]))*sin(spin1[2]) + (cnew-cold)*cos(spin1[2])  +  (snew*cos(spin0new[1]-spin2[1]) - sold*cos(spin0old[1]-spin2[1]))*sin(spin2[2]) + (cnew-cold)*cos(spin2[2])  +  (snew*cos(spin0new[1]-spin3[1]) - sold*cos(spin0old[1]-spin3[1]))*sin(spin3[2]) + (cnew-cold)*cos(spin3[2])  +  (snew*cos(spin0new[1]-spin4[1]) - sold*cos(spin0old[1]-spin4[1]))*sin(spin4[2]) + (cnew-cold)*cos(spin4[2])   +   d*(cold^2-cnew^2)
end

function EnergyDiffTheta(spin0old, theta, spin1, spin2, spin3, spin4, d) # calculate once sin and cos of old and new spin instead of 4, and once spin and cos of spin2-3-4-5 instead if 2
    sold, snew = sin(spin0old[2]), sin(theta)
    cold, cnew = cos(spin0old[2]), cos(theta)
    phi = spin0old[1]
    return (snew-sold)*(cos(phi-spin1[1])*sin(spin1[2]) + cos(phi-spin2[1])*sin(spin2[2]) + cos(phi-spin3[1])*sin(spin3[2]) + cos(phi-spin4[1])*sin(spin4[2])) + (cnew-cold)*(cos(spin1[2])+cos(spin2[2])+cos(spin3[2])+cos(spin4[2]))   +   d*(cold^2-cnew^2)
end

function EnergyDiffPhi(spin0old, phi, spin1, spin2, spin3, spin4) # calculate once sin and cos of old and new spin instead of 4, and once spin and cos of spin2-3-4-5 instead if 2
    phiold = spin0old[1]
    return sin(spin0old[2])*((cos(phi-spin1[1])-cos(phiold-spin1[1]))*sin(spin1[2])  +  (cos(phi-spin2[1])-cos(phiold-spin2[1]))*sin(spin2[2])  +  (cos(phi-spin3[1])-cos(phiold-spin3[1]))*sin(spin3[2])  +  (cos(phi-spin4[1])-cos(phiold-spin4[1]))*sin(spin4[2]))
end

function SingleFlip(Lattice::Array{Float32, 3}, i::Int, j::Int, T::Float32, L::Int, σ::Float32, pi32::Float32, PBC::Bool, d::Float32, p::Float32)      # Propose a new configuration and either accept or reject it
    newspin, WhichFlip = NewSpin(Lattice[:,i,j], σ, pi32, p)
    if PBC==false
        Δ = 0f0
        if i!=1
            Δ += Correlation(newspin,Lattice[:,i-1,j]) - Correlation(Lattice[:,i,j],Lattice[:,i-1,j]) # This is the right sign as the bounding energy is negative (and is calculate positively)
        end
        if i != L
            Δ += Correlation(newspin,Lattice[:,i+1,j]) - Correlation(Lattice[:,i,j],Lattice[:,i+1,j])
        end
        if j != 1
            Δ += Correlation(newspin,Lattice[:,i,j-1]) - Correlation(Lattice[:,i,j],Lattice[:,i,j-1])
        end
        if j != L
            Δ += Correlation(newspin,Lattice[:,i,j+1]) - Correlation(Lattice[:,i,j],Lattice[:,i,j+1])
        end
    else
        Δ = EnergyDiff(Lattice[:,i,j], newspin, Lattice[:,mod(i-2,L)+1,j],Lattice[:,i,mod(j-2,L)+1],Lattice[:,mod(i,L)+1,j],Lattice[:,i,mod(j,L)+1],d)
        # Δ = Correlation(newspin,Lattice[:,mod(i-2,L)+1,j]) - Correlation(Lattice[:,i,j],Lattice[:,mod(i-2,L)+1,j])   +   Correlation(newspin,Lattice[:,mod(i,L)+1,j]) - Correlation(Lattice[:,i,j],Lattice[:,mod(i,L)+1,j])   +   Correlation(newspin,Lattice[:,i,mod(j-2,L)+1]) - Correlation(Lattice[:,i,j],Lattice[:,i,mod(j-2,L)+1])   +   Correlation(newspin,Lattice[:,i,mod(j,L)+1]) - Correlation(Lattice[:,i,j],Lattice[:,i,mod(j,L)+1])
    end
    Δ += d * (cos(Lattice[2,i,j])^2 - cos(newspin[2])^2)
    if exp(Δ/T) > rand(TaskLocalRNG())
        Lattice[:,i,j] .= [mod(newspin[1], 2*pi32), newspin[2]]
        acceptance=true
        ΔE=-Δ
    else
        acceptance=false
        ΔE=0f0
    end
    return Lattice[:,i,j], ΔE, acceptance, WhichFlip
end

function NewSpin(spin::Vector, σ::Float32, pi32::Float32, p::Float32) # Heisenberg -- propose a new spin
    Whichflip = rand(TaskLocalRNG()) < p
    if Whichflip == true
        spin[2] = pi32 - spin[2]
    else
        alpha   = rand(TaskLocalRNG(), Float32)*2*pi32
        beta    = min(abs(rand(TaskLocalRNG(), Normal(0f0, σ))), pi32)
        spin[1] = atan(sin(spin[1])*(cos(spin[2])*sin(beta)*cos(alpha) + sin(spin[2])*cos(beta)) + cos(spin[1])*sin(beta)*sin(alpha),   cos(spin[1])*(cos(spin[2])*sin(beta)*cos(alpha) + sin(spin[2])*cos(beta)) - sin(spin[1])*sin(beta)*sin(alpha)) # alpha c'est le phi'
        spin[2] = acos(max(-1,min(1,cos(spin[2])*cos(beta) - sin(spin[2])*sin(beta)*cos(alpha))))
    end
    return spin, Whichflip+1
end

function SingleFlipXY(Lattice::Array{Float32, 3}, i::Int, j::Int, T::Float32, L::Int, σ::Float32, pi32::Float32, PBC::Bool, d::Float32, p::Float32)      # Propose a new configuration and either accept or reject it
    if rand(TaskLocalRNG()) > p
        newspin = [NewPhi(Lattice[1,i,j], σ[1], pi32), Lattice[2,i,j]]
        WhichAngle=1
    else
        newspin = [Lattice[1,i,j], NewTheta(Lattice[2,i,j], σ[2], pi32)]
        WhichAngle=2
    end
    
    Δ = 0f0
    if PBC==false
        if i!=1
            Δ += Correlation(newspin,Lattice[:,i-1,j]) - Correlation(Lattice[:,i,j],Lattice[:,i-1,j]) # This is the right sign as the bounding energy is negative (and is calculate positively)
        end
        if i != L
            Δ += Correlation(newspin,Lattice[:,i+1,j]) - Correlation(Lattice[:,i,j],Lattice[:,i+1,j])
        end
        if j != 1
            Δ += Correlation(newspin,Lattice[:,i,j-1]) - Correlation(Lattice[:,i,j],Lattice[:,i,j-1])
        end
        if j != L
            Δ += Correlation(newspin,Lattice[:,i,j+1]) - Correlation(Lattice[:,i,j],Lattice[:,i,j+1])
        end
    else
        Δ += Correlation(newspin,Lattice[:,mod(i-2,L)+1,j]) - Correlation(Lattice[:,i,j],Lattice[:,mod(i-2,L)+1,j])   +   Correlation(newspin,Lattice[:,mod(i,L)+1,j]) - Correlation(Lattice[:,i,j],Lattice[:,mod(i,L)+1,j])   +   Correlation(newspin,Lattice[:,i,mod(j-2,L)+1]) - Correlation(Lattice[:,i,j],Lattice[:,i,mod(j-2,L)+1])   +   Correlation(newspin,Lattice[:,i,mod(j,L)+1]) - Correlation(Lattice[:,i,j],Lattice[:,i,mod(j,L)+1])
    end
    Δ += d * (cos(Lattice[2,i,j])^2 - cos(newspin[2])^2)
    if exp(Δ/T) > rand(TaskLocalRNG())
        Lattice[:,i,j] .= [mod(newspin[1], 2*pi32), newspin[2]]
        acceptance=true
        ΔE=-Δ
    else
        acceptance=false
        ΔE=0f0
    end
    return Lattice[:,i,j], ΔE, acceptance, WhichAngle
end

function NewPhi(phi, σϕ, pi32)
    return mod(rand(TaskLocalRNG(), Normal(phi, σϕ)), 2*pi32)
end

function NewTheta(theta, σθ, pi32)
    theta = mod(rand(TaskLocalRNG(), Normal(theta, σθ)), 2*pi32)
    if theta > pi32; theta = 2*pi32 - theta
    elseif theta < 0; theta = - theta
    end
    return  theta
end

function LatticeSweep(Lattice::Array{Float32, 3}, L::Int, E::Float32, T::Float32, σ::Float32, pi32::Float32, PBC::Bool, d::Float32, p::Float32)    #MH lattice sweep
    Accept = [0,0]
    Tr = [0,0]
    for j=1:L
        for k=1:L
            Lattice[:,j,k], ΔE, a, Whichflip = SingleFlip(Lattice,j,k,T,L,σ,pi32, PBC, d, p)
            Accept[Whichflip]+=a
            Tr[Whichflip]+=1
            E+=ΔE
        end
    end
    return Lattice, E, Accept, Tr
end

function MultipleIsingFlips(Lattice::Array{Float32, 3}, L::Int, T::Float32, pi32::Float32, d::Float32, PBC::Bool)      # Wolff update
    p = 1 - ℯ^(-2/T)    # what should we do about d ?
    k, l = rand(TaskLocalRNG(), 1:L, 2)
    InCluster = falses(L,L)
    InCluster[k,l] = true
    CheckNeighbor = [[k,l]]
    while CheckNeighbor != []
        i,j = pop!(CheckNeighbor)
        NewCheckNeighbor, InCluster = AddNeighbors(i,j, L, InCluster, p, Lattice[2,:,:],pi32) # check first i-1,j; then i,j-1; i+1,j and finally i,j+1
        CheckNeighbor = unique(vcat(CheckNeighbor,NewCheckNeighbor))
    end
    Lattice[2,:,:] = InCluster .* (pi32 .- Lattice[2,:,:]) .+ .!InCluster .* Lattice[2,:,:]
    E = Energy(Lattice, L, PBC, d)
    return Lattice, E
end

function AddNeighbors(i::Int,j::Int, L::Int, InCluster::BitMatrix, p::Float32, ZLattice::Array{Float32, 2}, pi32::Float32)    # add or not to the cluster the neighbors of a site already in the cluster
    NewCheckNeighbor = []
    Zsign = sign(ZLattice[i,j]-pi32/2)
    neighbor = [[mod(i-2,L)+1,j],[i,mod(j-2,L)+1],[mod(i,L)+1,j],[i,mod(j,L)+1]]    # The four neighbors
    for n=1:4                                       # To check the four neighbors
        k,l = neighbor[n]
        if InCluster[k,l] == false && sign(ZLattice[k,l]-pi32/2) == Zsign && rand(TaskLocalRNG()) < p
            InCluster[k,l] = true
            push!(NewCheckNeighbor, [k,l])
        end
    end
    return NewCheckNeighbor, InCluster
end

function Swap(Replicas::Vector{Replica}, i::Int)    # Try to swap with the next temperature
    j=i+1
    if (1/Replicas[i].T - 1/Replicas[j].T)*(Replicas[i].E - Replicas[j].E) > rand(TaskLocalRNG())
        Replicas[i].Lattice, Replicas[j].Lattice = Replicas[j].Lattice, Replicas[i].Lattice
        Replicas[i].E,       Replicas[j].E       = Replicas[j].E,       Replicas[i].E
        Replicas[i].σ,       Replicas[j].σ       = Replicas[j].σ,       Replicas[i].σ
        return true
    else
        return false
    end
end

function MHStep(step::Int, rep::Replica, L::Int, d::Float32, p::Float32, PBC::Bool, pi32::Float32)
    IsingFlip = mean(cos.(rep.Lattice[2,:,:]).^2) > .75     # require that in average [theta < pi/6 or theta > 5 pi/6] to have the possibility (with proba p) to Ising-flip
    if mod(step,1000)==0 && IsingFlip==true
        rep.Lattice, rep.E = MultipleIsingFlips(rep.Lattice, L, rep.T, pi32, d, PBC)    # Wolff algorithm
    else
        rep.Lattice, rep.E, Accept, Tr = LatticeSweep(rep.Lattice, L, rep.E, rep.T, rep.σ, pi32, PBC, d, p*IsingFlip)
        if Tr[1] > 20
            rep.σ = min(max(rep.σ*2*Accept[1]/Tr[1], 1/sqrt(step+1000)), 10)     # if the acceptance is higher than .5, sigma should increase (decreased) to explore more (less) the phasespace to get closer to .5
            rep.Acceptance += Accept
            rep.Try += Tr
        end
    end
    return rep
end

function MH_parallel_tempering(L::Int, N::Int, T::Vector{Float32}, burn::Int, pi32::Float32, AllLattices::Vector{Int}, d::Float32, p::Float32, Save::Bool=false, Skip::Int=10, swap::Int=101)     # Sampler
    nT = length(T)
    Replicas = Vector{Replica}(undef, nT)
    for i=1:nT  # Initialize Replicas
        lat = Initial_lattice(L, pi32)
        Replicas[i] = Replica(lat, Energy(lat, L, PBC, d), T[i], .25f0, [0,0], [0,0])
    end
    Nswap = round(Int, N/swap)
    Nmeasurement = round(Int, (N - burn)/Skip)
    energies = zeros(Float32, nT, Nmeasurement)
    # Energies2 = zeros(Float32, nT, Int(Nmeasurement/Skip))
    mag  = zeros(Float32, nT, Nmeasurement)
    correl = zeros(Float32, nT, length(RowCorr(Replicas[1].Lattice, L, PBC)))
    SwapAcceptance = zeros(Int, nT-1)
    
    for i=1:Nswap
        Threads.@threads for n=1:nT
            for j=1:swap
                k = (i-1)*swap+j
                Replicas[n] = MHStep(k, Replicas[n], L, d, p, PBC, pi32)
                if mod(k, Skip) == 0 && k > burn            # to not measure every lattice sweeps
                    lat = Replicas[n].Lattice
                    m = round(Int,(k-burn)/Skip)
                    energies[n, m] = copy(Replicas[n].E)
                    mag[n, m] = sqrt(mean(cos(phi)*sin(theta) for phi in lat[1,:,:], theta in lat[2,:,:])^2 + mean(sin(phi)*sin(theta) for phi in lat[1,:,:], theta in lat[2,:,:])^2 + mean(cos(theta) for theta in lat[2,:,:])^2)
                    correl[n,:] += RowCorr(lat, L, PBC)
                    # Energies2[n, m] = copy(Energy(lat, L, PBC, d))/(!PBC*L*(L-1) + PBC*L^2)
                end
            end
        end
        for n=1:nT-1
            SwapAcceptance[n] += Swap(Replicas, n)
        end
    end
    SwapAccept = SwapAcceptance/Nswap
    energies /= !PBC*L*(L-1) + PBC*L^2
    if Save==true
        for n=1:nT
            Energies, Mag, corr = energies[n,:], mag[n,:], correl[n] 
            accept = Replicas[n].Acceptance ./ Replicas[n].Try
            @save "Data/$(AllLattices)_$N/$(L)_$(T[n])_$d.jld2" Energies Mag corr accept# Energies2
        end
        @save "Data/$(AllLattices)_$N/swap_$(L)_$d.jld2" SwapAccept
    end
    return mean(energies[1,:]), #=mean(Energies2[1,:]),=# mean(mag[1,:])
end

function MHvideo(Lattice::Array{Float32, 3}, L::Int, N::Int, T::Float32, Nbin::Int, burn::Int, pi32::Float32, PBC::Bool, d::Float32, p::Float32)     # Sampler
    T0 = 1.3f0
    σ=.25
    acceptance=[0,0]                  # around, Ising, multiflips
    Try = [0,0]
    Energies= zeros(Float32, N-burn)
    Mag  = zeros(Float32, N-burn)
    E = Energy(Lattice, L, PBC, d)
    corr = zeros(Float32, length(RowCorr(Lattice, L, PBC)))
    Lattices=[]
    
    for i=1:N                   # Calculate the energy at each flip
        Accept = [0,0]
        Tr = [0,0]
        IsingFlip = mean(cos.(Lattice[2,:,:]).^2) > .75     # require that in average [theta < pi/6 or theta > 5 pi/6] to have the possibility (with proba p) to Ising-flip
        if T0 > T
            T0 = 1.3f0*(1f0 - Float32(i/burn))
        else; T0 = T
        end
        if mod(i,40)==0 && IsingFlip==true                   # Wolff algorithm
            Lattice, E = MultipleIsingFlips(Lattice, L, T, pi32, d, PBC)
        else
            Lattice, E, Accept, Tr = LatticeSweep(Lattice, L, E, Accept, Tr, T0, σ, pi32, PBC, d, p*IsingFlip)
            acceptance += Accept
            Try += Tr
            if Tr[1] > 3
                σ = min(max(σ*2*Accept[1]/Tr[1], 1/sqrt(i+1000)), 2)     # if the acceptance is higher than .5, sigma should increase else to be decreased to explore more the phasespace
            end
        end
        if mod(i,Int(N/200))==0
            push!(Lattices, copy(Lattice))
        end
        if i > burn
            j = i-burn
            Energies[j] = E
            Mag[j] = sqrt(mean(cos(phi)*sin(theta) for phi in Lattice[1,:,:], theta in Lattice[2,:,:])^2 + mean(sin(phi)*sin(theta) for phi in Lattice[1,:,:], theta in Lattice[2,:,:])^2 + mean(cos(theta) for theta in Lattice[2,:,:])^2)
            corr += RowCorr(Lattice, L, PBC)
        end
    end
    println("end of calculation phase")
    Energies /= !PBC*L*(L-1) + PBC*L^2
    Magbin = Binor(Mag, Nbin)
    # χmean     = (mean(Mag.^2)-mean(Mag)^2)/T*L^2
    # Cmean     = (mean(Energies.^2)-mean(Energies)^2)/T^2*L^2
    EnergyBins=Binor(Energies, Nbin)
    println(typeof(Lattices))

    anim = @animate for j in eachindex(Lattices)
        heatmap(matrixcolor(Lattices[j], L), aspect_ratio = 1, size = (400,400), colormap = :coolwarm, legend = false, framestyle=:box, title = "T=$T  d=$d  $(Int(j*N/200))")    
    end
    # anim = @animate for j in eachindex(Lattices)
    #     p1 = heatmap(matrixcolor(Lattices[j], L),
    #         aspect_ratio = 1, size = (288,288),
    #         colormap = :coolwarm, legend = false,
    #         framestyle=:box, title = "matrixcolor")
    #     p2 = heatmap(matrixcolorXY(Lattices[j], L),
    #         aspect_ratio = 1, size = (288,288),
    #         colormap = :coolwarm, legend = false,
    #         framestyle=:box, title = "matrixcolorXY")
    #     p3 = heatmap(matrixcolorZ(Lattices[j], L),
    #         aspect_ratio = 1, size = (288,288),
    #         colormap = :coolwarm, legend = false,
    #         framestyle=:box, title = "matrixcolorZ")
    #     plot(p1, p2, p3, layout = (1,3), size = (864,288), title = "T=$T   $j")
    # end

    gif(anim, "C:/Users/yanis/Desktop/Cours/These/Heisenberg/Spin_$N.gif", fps=10)

    return mean(Energies), std(EnergyBins), mean(Mag), std(Magbin)#=, χmean, Errorpropagation(Magbin,std(Magbin))/T*L^2 #=std(χbin)=#, Cmean, Errorpropagation(EnergyBins,std(EnergyBins))/T^2*L^2, corr/(N-burn+1)=#, acceptance ./ Try
end

function MH(Lattice::Array{Float32, 3}, L::Int, N::Int, T::Float32, burn::Int, pi32::Float32, AllLattices::Vector, d::Float32, p::Float32, Save::Bool, Skip=10::Int)     # Sampler
    T0 = 1.3f0
    σ=.25f0
    Nmeasurement = N - burn
    acceptance = [0,0]              # around, Ising, multiflips
    Try = [0,0]
    Accept = [0,0]                  # accept is usefull to keep track of recent the acceptance rate whereas acceptance keep track of the total
    Tr = [0,0]                      # same reason as accept
    E = Energy(Lattice, L, PBC, d)
    Energies = zeros(Float32, Int(Nmeasurement/Skip))
    # Energies2 = zeros(Float32, Int(Nmeasurement/Skip))
    Mag  = zeros(Float32, Int(Nmeasurement/Skip))
    corr = zeros(length(RowCorr(Lattice, L, PBC)))

    for i=1:N                   # Calculate the energy at each flip
        IsingFlip = mean(cos.(Lattice[2,:,:]).^2) > .75     # require that in average [theta < pi/6 or theta > 5 pi/6] to have the possibility (with proba p) to Ising-flip
        if T0 > T                                           # To begin with a "high" temperature to help it converges more easily, and decrease the temperature until the wanted one
            T0 = 2.8f0*(1 - i/burn)
        else; T0 = T
        end
        if mod(i,1000)==0 && IsingFlip==true                # Wolff algorithm
            Lattice, E = MultipleIsingFlips(Lattice, L, T, pi32, d, PBC)
        else
            Lattice, E, Accept, Tr = LatticeSweep(Lattice, L, E, Accept, Tr, T0, σ, pi32, PBC, d, p*IsingFlip)
            if Tr[1] > 1000
                σ = min(max(σ*2*Accept[1]/Tr[1], 1/sqrt(i+1000)), 10)     # if the acceptance is higher than .5, sigma should increase (decreased) to explore more (less) the phasespace to get closer to .5
                acceptance += Accept
                Try += Tr
                Accept = [0,0]
                Tr = [0,0]
            end
        end
        if mod(i, Skip) == 0 && i > burn            # to not measure every lattice sweeps
            j = Int((i-burn)/Skip)
            Energies[j] = copy(E)
            Mag[j] = sqrt(mean(cos(phi)*sin(theta) for phi in Lattice[1,:,:], theta in Lattice[2,:,:])^2 + mean(sin(phi)*sin(theta) for phi in Lattice[1,:,:], theta in Lattice[2,:,:])^2 + mean(cos(theta) for theta in Lattice[2,:,:])^2)
            corr += RowCorr(Lattice, L, PBC)
            # Energies2[j] = copy(Energy(Lattice, L, PBC, d))/(!PBC*L*(L-1) + PBC*L^2)
        end
    end
    Energies /= !PBC*L*(L-1) + PBC*L^2
    accept = acceptance ./ Try
    if Save==true; @save "Data/$(AllLattices)_$N/$(L)_$(T)_$d.jld2" Energies Mag corr accept# Energies2
    end
    return mean(Energies), #=mean(Energies2),=# mean(Mag)
end

function MHXY(Lattice, L, N, T, burn, pi32, AllLattices, d, Save, Skip=10)     # Sampler
    T0 = 1.3f0
    σ=[.8f0, .45f0]
    p=.5
    Nmeasurement = N - burn
    acceptance = [0,0]              # around, Ising, multiflips
    Try = [0,0]
    E = Energy(Lattice, L, PBC, d)
    Energies = zeros(Float32, Int(Nmeasurement/Skip))
    Mag  = zeros(Float32, Int(Nmeasurement/Skip))
    corr = zeros(length(RowCorr(Lattice, L, PBC)))
    magz, absmagz = zeros(Float32, Int(Nmeasurement/Skip)), zeros(Float32, Int(Nmeasurement/Skip))
    sigma1, sigma2 = zeros(Float32, Int(Nmeasurement/Skip)), zeros(Float32, Int(Nmeasurement/Skip))

    for i=1:N                   # Calculate the energy at each flip
        Accept = [0,0]
        Tr = [0,0]
        if T0 > T                                           # To begin with a "high" temperature to help it converges more easily, and decrease the temperature until the wanted one
            T0 = 1.3f0*(1f0 - Float32(round(i/N;digits=1)))
        else; T0 = T
        end
        Lattice, E, Accept, Tr = LatticeSweep(Lattice, L, E, Accept, Tr, T0, σ, pi32, PBC, d, 0)
        acceptance += Accept
        Try += Tr
        σ = min.(σ * 2 .*Accept./Tr .+.0000001f0, [2,1])     # if the acceptance is higher than .5, sigma should increase else to be decreased to explore more the phasespace
        p = σ[2]/sum(σ)                             # if σtheta (σ[2]) is very low (and M_z around 0), it means that we are probably in XY limit, thus it will be more intersting to move around phi as theta already stable.
        if mod(i, Skip) == 0 && i > burn            # to not measure every lattice sweeps
            j = Int((i-burn)/Skip)
            Energies[j] = copy(E)
            Mag[j] = sqrt(mean(cos(phi)*sin(theta) for phi in Lattice[1,:,:], theta in Lattice[2,:,:])^2 + mean(sin(phi)*sin(theta) for phi in Lattice[1,:,:], theta in Lattice[2,:,:])^2 + mean(cos(theta) for theta in Lattice[2,:,:])^2)
            corr += RowCorr(Lattice, L, PBC)
            magz[j] = mean(cos(theta) for theta in Lattice[2,:,:])
            sigma1[j] = σ[1]
            sigma2[j] = σ[2]
            absmagz[j] = mean(abs(cos(theta)) for theta in Lattice[2,:,:])
        end
    end
    Energies /= !PBC*L*(L-1) + PBC*L^2
    accept = acceptance ./ Try
    if Save==true; @save "Data/$(AllLattices)_$N/$(L)_$(T)_$d.jld2" Energies Mag corr accept
    end
    return Energies[end], magz, accept, σ, sigma1, sigma2, absmagz
end

function MHCorr(Lattice, L, N, T, burn, pi32, AllLattices, d, p, Save, Skip=10)     # Sampler
    T0 = 1.3f0
    σ=.25f0
    Nmeasurement = N - burn
    acceptance = [0,0]              # around, Ising, multiflips
    Try = [0,0]
    Accept = [0,0]                  # accept is usefull to keep track of recent the acceptance rate whereas acceptance keep track of the total
    Tr = [0,0]                      # same reason as accept
    E = Energy(Lattice, L, PBC, d)
    Energies = zeros(Float32, Int(Nmeasurement/Skip))
    Energies2 = zeros(Float32, Int(Nmeasurement/Skip))
    Mag  = zeros(Float32, Int(Nmeasurement/Skip))
    Lattices=[]
    corr = zeros(length(RowCorr(Lattice, L, PBC)))
    
    for i=1:N                   # Calculate the energy at each flip
        IsingFlip = mean(cos.(Lattice[2,:,:]).^2) > .75     # require that in average [theta < pi/6 or theta > 5 pi/6] to have the possibility (with proba p) to Ising-flip
        if T0 > T                                           # To begin with a "high" temperature to help it converges more easily, and decrease the temperature until the wanted one
            T0 = 2.8f0*(1 - Float32(round(i/burn;digits=1)))
        else; T0 = T
        end
        if mod(i,1000)==0 && IsingFlip==true                # Wolff algorithm
            Lattice, E = MultipleIsingFlips(Lattice, L, T, pi32, d, PBC)
        else
            Lattice, E, Accept, Tr = LatticeSweep(Lattice, L, E, Accept, Tr, T0, σ, pi32, PBC, d, p*IsingFlip)
            if Tr[1] > 1000
                σ = min(max(σ*2*Accept[1]/Tr[1], 1/sqrt(i+1000)), 10)     # if the acceptance is higher than .5, sigma should increase (decreased) to explore more (less) the phasespace to get closer to .5
                acceptance += Accept
                Try += Tr
                Accept = [0,0]
                Tr = [0,0]
            end
        end
        if mod(i, Skip) == 0 && i > burn            # to not measure every lattice sweeps
            j = Int((i-burn)/Skip)
            Energies[j] = copy(E)
            Mag[j] = sqrt(mean(cos(phi)*sin(theta) for phi in Lattice[1,:,:], theta in Lattice[2,:,:])^2 + mean(sin(phi)*sin(theta) for phi in Lattice[1,:,:], theta in Lattice[2,:,:])^2 + mean(cos(theta) for theta in Lattice[2,:,:])^2)
            corr += RowCorr(Lattice, L, PBC)
            Energies2[j] = copy(Energy(Lattice, L, PBC, d))
            push!(Lattices, copy(Lattice))
        end
    end
    Energies /= !PBC*L*(L-1) + PBC*L^2
    Energies2 /= !PBC*L*(L-1) + PBC*L^2
    accept = acceptance ./ Try
    if Save==true; @save "Data/$(AllLattices)_$N/$(L)_$(T)_$d.jld2" Energies Mag corr accept Energies2 Lattices
    end
    return mean(Energies), mean(Energies2), mean(Mag)
end

function matrixcolor(Lattice::Array{Float32, 3}, L::Int)     # represent the phase by a color
    B=zeros(RGB{Float32},L,L)
    for i=1:L
        for j=1:L
            B[i,j]= RGB{Float32}(.5+cos(Lattice[1,i,j])*sin(Lattice[2,i,j])/2, .5+cos(Lattice[2,i,j])/2, .5+sin(Lattice[1,i,j])*sin(Lattice[2,i,j])/2)
        end
    end
    return B
end

function matrixcolorXY(Lattice, L)     # represent the phase by a color
    B=zeros(RGB{Float32},L,L)
    for i=1:L
        for j=1:L
            B[i,j]= RGB{Float32}(.5+cos(Lattice[1,i,j])/2, 0, .5+sin(Lattice[1,i,j])/2)
        end
    end
    return B
end

function matrixcolorZ(Lattice, L)     # represent the phase by a color
    B=zeros(RGB{Float32},L,L)
    for i=1:L
        for j=1:L
            B[i,j]= RGB{Float32}(0, .5+cos(Lattice[2,i,j])/2, 0)
        end
    end
    return B
end

function RowCorr(Lattice::Array{Float32, 3}, L::Int, PBC::Bool)   # Return the average row Correlation of the matrix (in a Correlation vector : first neigbour, second neigbour...)
    corr=[]
    if PBC==false
        for i=1:(L-1)
            a=0
            for j=1:L
                for k=1:(L-i)
                    a+=Correlation(Lattice[:,j,k],Lattice[:,j,mod(k+i-1,L)+1])
                end
            end
            push!(corr, a/(L-i)/L)
        end
    else
        for i=1:round(Int, L/2-1)
            a=0
            for j=1:L
                for k=1:L
                    a+=Correlation(Lattice[:,j,k],Lattice[:,j,mod(k+i-1,L)+1]) #cos(Lattice[1,j,k]-Lattice[1,j,mod(k+i-1,L)+1]) * cos(Lattice[2,j,k]-Lattice[2,j,mod(k+i-1,L)+1])
                end
            end
            push!(corr, a)
        end
    end
    return corr/L^2
end


#------     Analysis


function Get_T_d(N,L)  # for the first lattice size (supposing they all have the same T and d), gives which T have each d
    files = glob("$(L[1])*.jld2", ("Data/$(L)_$N"))
    for f in files
        name = splitext(last(splitpath(f)))[1]
        parts = split(name, "_")
        if length(parts) != 3
            @warn "Bad filename format" f parts
        end
    end

    parsed = map(files) do f
        name = splitext(last(splitpath(f)))[1]   # robust path handling
        nums = parse.(Float32, split(name, "_"))
        (; n2 = nums[2], n3 = nums[3])
    end

    # all_L = [p.n1 for p in parsed]  # Vector of all L
    all_d = unique(sort([p.n3 for p in parsed]))    # Vector of all d
    all_T = unique(sort([p.n2 for p in parsed]))    # Vector of all T
    T_for_d = Dict{Float32, Vector{Float32}}() # Dictionary mapping each #3 → list of #2
    for p in parsed
        push!(get!(T_for_d, p.n3, Float32[]), p.n2)
    end
    T_for_d = Dict(k => sort(v) for (k,v) in T_for_d)
    return all_T, all_d, T_for_d
end

function Binor(arr, Nbin) # take them Nbin by Nbin and average each the bins.
    Bins=[]
    Nperbin = round(Int, length(arr)/Nbin)
    for i=1:Nbin
        push!(Bins, mean(arr[((i-1)*Nperbin)+1:i*Nperbin]))
    end
    return Bins
end

function Binor2nd(arr, Nbin, T, L) # USELESS FOR MH (as I use error propag)
    Bins=[]
    Nperbin = round(Int, length(arr)/Nbin)
    for i=1:Nbin
        push!(Bins, (mean(arr[(i-1)*Nperbin+1:i*Nperbin].^2) - mean(arr[(i-1)*Nperbin+1:i*Nperbin])^2)/T*L^2)
    end
    return Bins
end

function basicplot1(L, T_for_d, dvect, d, vect, ytitle="", error=[], save=false)   # read data from a file and plot it
    i = findfirst(==(d), dvect)
    endT = length(vect[1,:,i])
    while endT >= 1 && vect[1,endT,i]==0
        endT -= 1
    end
    p=plot()
    if typeof(vect) == Matrix{Any}
        for t=1:max(round(Int, length(T_for_d[d])/8),1):length(T_for_d[d])
            plot!(collect(1:Int(L[end]/2-1)), vect[end,:][t], label="$(T[t])",legend=:topright) # yet not update to work
        end
        title!(p, "Correlation (L=$(L[end])) as a function of distance")
        xlabel!(p, "Distance (site)")
    else
        if error != []
            if ytitle=="C" || ytitle=="Susc"; lim=(max(0, minimum(vect[:,:,i])),maximum(vect[:,3:end,i])+maximum(error[:,3:end,i]))
            else; lim = :auto
            end
            for l in eachindex(L)
                plot!(T_for_d[d], vect[l, 1:endT,i], yerr = Vector(error[l,1:endT,i]), markerstrokecolor=:auto, label="$(L[l])", ylims = lim)
            end
        else
            for l in eachindex(L)
                plot!(T_for_d[d], vect[l, 1:endT,i], label="$(L[l])", ylims=(max(0, minimum(vect[:,1:endT,i])),maximum(vect[:,1:endT,i])))
            end
        end
        title!(p, ytitle*" as a function of T with d = $d")
        xlabel!(p, "Temperature")
    end
    ylabel!(p, ytitle)
    if save ==true
        savefig("Plot/"*ytitle*".pdf")
    end
    println("size vect: ", size(vect))
    println("size vect[end,:,i]: ", size(vect[end,:,i]))
    println("endT: ", endT)
    display(p)
end

function basicplotd(T_for_d, d, vect, ytitle="", error=[], i=length(vect[:,1,1]), save=false)   # read data from a file and plot it
    p=plot()
    pal = cgrad([RGB(.4,.6,1), RGB(0,0,.5), RGB(0,0,0), RGB(.6,0,0), RGB(1,.6,.6)], [0., .4999, .5, .5001, 1.], categorical = false)
    if error != [] 
        if ytitle=="C"#= || ytitle=="Susc"=#; lim=(max(0, minimum(vect[i,:,:])),maximum(vect[i,3:end,:])+maximum(error[i,3:end,:]))
        else; lim = :auto
        end
        for z in eachindex(d)
            rescaleE=0          # To rescale the Energy (only when y = Energy)
            endT = length(vect[i,:,z])
            while endT >= 1 && vect[i,endT,z]==0
                endT -= 1
            end
            if ytitle=="E" && d[z]<0; rescaleE = d[z]; end
            x = PlotColord(d[z], d[1], d[end])
            plot!(T_for_d[d[z]], vect[i, 1:endT, z].-rescaleE, yerr = Vector(error[i, 1:endT, z]), label="$(d[z])", ylims=lim, seriescolor = pal[x], linecolor = pal[x], markercolor = pal[x], markerstrokecolor = pal[x], ecolor = pal[x])
        end
    else
        for z in eachindex(d)
            endT = length(vect[i,:,z])
            while endT >= 1 && vect[i,endT,z]==0
                endT -= 1
            end
            plot!(T_for_d[d[z]], vect[i, 1:endT, z], label="$(d[z])", seriescolor = pal[z])
        end
    end
    title!(p, ytitle*" as a function of T")
    xlabel!(p, "Temperature")
    ylabel!(p, ytitle)
    if save ==true
        savefig("Plot/"*ytitle*".pdf")
    end
    display(p)
end

function PlotColord(x, dmin, dmax)
    if x < 0
        return .5 - .5* sqrt(x / dmin)
    elseif x > 0
        return .5 + .5 * sqrt(x / dmax)
    else
        return .5
    end
end

function Errorpropagation(v, Δ) # propagation of error for c and χ from the error on E and M
    return sqrt(sum((v .- mean(v)).^2)/length(v))*2*Δ
end

function interpmax(L, T, y) # interpolate and find the max
    Max = []
    xinterp = collect(.1:.001:2)
    p=plot()
    for l in eachindex(L)
        yinterp=Spline1D(T, y[l,:], k=3)(xinterp)
        push!(Max, (findmax(yinterp[301:end])[1], xinterp[findmax(yinterp[301:end])[2]]))
    end
    return Max
end

function linear(x, p) # linear function
    return p[1]*x + p[2]   
end

function crit(L, T, y, graph="") # calculate α & γ
    ymax = interpmax(L, T, y)
    println(ymax)
    fit =  curve_fit(linear, log.(L), log.(first.(ymax)), if (y=="c"); [.1,-.7]; elseif (y=="susc"); [2.,-4.]; else [.9,-.9] end)
    println(typeof(fit), "\n", fit)
    m,p = coef(fit)
    σ = stderror(fit)
    if graph!=""
        a=plot()
        plot!(a, log.(L), log.(first.(ymax)), seriestype=:scatter)
        plot!(a, log.(L), m*log.(L) .+p)
        xlabel!(a, "Ln(Lattice length)")
        ylabel!(a, "Ln(max("*graph*"))")
        display(a)
    end
    return m, σ[1]
end

function critlength(T, data, t, ln=false) # Calculate critical exp or correlation length
    closestT = argmin(abs.(T .- t)); println(T[closestT])
    data=data[Int(length(Corrmean)/length(T))*closestT]
    negativ = findfirst(x -> x < 1e-5, data)
    if typeof(negativ) == Nothing; l = length(data); else; l=negativ-1; end
    println(l)
    x=collect(1:l)
    fit =  curve_fit(linear, log.(x)*ln + !ln*x, log.(data[1:l]), [-.09,9.])
    m,p = coef(fit); σ = stderror(fit)

    a=plot()
    plot!(a, log.(x)*ln + !ln*x, log.(data[1:l]), seriestype=:scatter)
    plot!(a, log.(x)*ln + !ln*x, m*(log.(x)*ln + !ln*x) .+p)
    ylabel!(a, "Ln(Correlation)")
    if ln == true; xlabel!(a, "Ln(Distance)"); ξ=σξ=0
    else; xlabel!(a, "Distance"); ξ =-1/m; σξ = σ[1]/m^2;end
    title!("Correlation vs distance at T = $(T[closestT])")
    display(a)
    return m*ln + !ln*ξ, σ[1]*ln+!ln*σξ
end

function MeanCorrTime(Lattices)
    L = length(Lattices[1][:,1,1])
    n = length(Lattices)
    Corr = zeros(n-1)
    for k=1:n-1
        for l=k+1:n
            C = 0
            for i=1:L
                for j=1:L
                    C += Correlation(Lattices[k][:,i,j],Lattices[l][:,i,j])
                end
            end
            Corr[l-k] += abs(C)/(n-l+k)
        end
    end
    return Corr*L^(-2)
end

function CorrPlot(AllLattices, T_for_d, dvect, d, L, l,save=false)
    p=plot()
    pal = cgrad([:lightblue, :green, :red])
    for i in eachindex(T_for_d[dvect[d]])
        MCT = MeanCorrTime(AllLattices[l,i,d])
        x = T_for_d[dvect[d]][i]/T_for_d[dvect[d]][end]
        plot!(collect(1:length(MCT))*10, MCT, label="$(T_for_d[dvect[d]][i])", seriescolor = pal[x], linecolor = pal[x], markercolor = pal[x], markerstrokecolor = pal[x], ecolor = pal[x])
    end
    title!(p, "Correlation time, L=$(L[l]), d=$(dvect[d])")
    if save ==true
        savefig("Plot/CorrelationTime$(dvect[d]).pdf")
    end
    display(p)
end

# α : specific heat                     Ising : 0       XY : NOP Essential singularity
# β : zero field mag                    Ising : 1/8     XY : NOP No magnetization (To have a nonzero 𝑀, the correlation function must approach a constant at large distance. But in the 2D XY model: For 𝑇>𝑇BKT: correlations decay exponentially. For 𝑇<𝑇BKT: correlations decay as a power law)
# γ : zero field isothermal suscxx      Ising : 7/4     XY : NOP Essential divergence
# δ : Critical isothermal               Ising : 15      XY : 15
# ν : corr length                       Ising : 1       XY : NOP 𝜉 diverges exponentially

# continuous symmetries cannot be spontaneously broken in 2D