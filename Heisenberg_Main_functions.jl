mutable struct Replica
    Lattice::Array{Float32,3}
    T::Float32                  # Constant
    E::Float32                  # Instantaneous
    σ::Float32                  # Instantaneous
    σϕ::Float32                 # Instantaneous
    σθ::Float32                 # Instantaneous
    Acceptance::Vector{Int64}   # History
    Try::Vector{Int64}          # History
end

function Energy(Lattice::Array{Float32, 3}, L::Int, PBC::Bool, d::Float32)    # Calculate the Energy of a configuration
    E=0f0
    n = PBC ? L : L-1
    @inbounds for i=1:n
        i2 = i==L ? 1 : i+1
        for j=1:n
            j2 = j==L ? 1 : j+1
            E += -Correlation(Lattice[1, i,  j], Lattice[2, i, j], Lattice[1, i, j2], Lattice[2, i, j2]) - Correlation(Lattice[1, i,  j], Lattice[2, i,  j], Lattice[1, i2,  j], Lattice[2, i2, j]) + d * cos(Lattice[2, i, j])^2
        end
    end
    return E
end #   @inbounds almost divide the execution time by 2

@inline function Correlation(ϕ1::Float32, θ1::Float32, ϕ2::Float32, θ2::Float32)
    sin(θ1)*sin(θ2)*cos(ϕ1-ϕ2) + cos(θ1)*cos(θ2)
end

function Initial_lattice(L::Int, pi32::Float32)
    lattice = Array{Float32}(undef, 2, L, L)
    lattice[1, :, :] .= 2*pi32*rand(TaskLocalRNG(), Float32, L, L)
    lattice[2, :, :] .= acos.(1 .- 2*rand(TaskLocalRNG(), Float32, L, L))
    return lattice
end

@inline function EnergyDiff(spin0old::Vector{Float32}, spin0new::Vector{Float32}, spin1::Vector{Float32}, spin2::Vector{Float32}, spin3::Vector{Float32}, spin4::Vector{Float32}, d::Float32) # calculate the bounding energy difference between a spin and its neigbors and another one with the same neighbors
    sold, snew = sin(spin0old[2]), sin(spin0new[2])
    cold, cnew = cos(spin0old[2]), cos(spin0new[2])
    return (snew*cos(spin0new[1]-spin1[1]) - sold*cos(spin0old[1]-spin1[1]))*sin(spin1[2]) + (cnew-cold)*cos(spin1[2])  +  (snew*cos(spin0new[1]-spin2[1]) - sold*cos(spin0old[1]-spin2[1]))*sin(spin2[2]) + (cnew-cold)*cos(spin2[2])  +  (snew*cos(spin0new[1]-spin3[1]) - sold*cos(spin0old[1]-spin3[1]))*sin(spin3[2]) + (cnew-cold)*cos(spin3[2])  +  (snew*cos(spin0new[1]-spin4[1]) - sold*cos(spin0old[1]-spin4[1]))*sin(spin4[2]) + (cnew-cold)*cos(spin4[2])   +   d*(cold^2-cnew^2) # old - new for the last term, as H = - Σ Si Sj + d Σ Sz^2 (sign difference between these two terms)
end

@inline function EnergyDiffTheta(spin0old::Vector{Float32}, theta::Float32, spin1::Vector{Float32}, spin2::Vector{Float32}, spin3::Vector{Float32}, spin4::Vector{Float32}, d::Float32) # Like EnergyDiff but optimized for a change in theta only (for XY-theta sweep or Ising flip)
    sold, snew = sin(spin0old[2]), sin(theta)
    cold, cnew = cos(spin0old[2]), cos(theta)
    phi = spin0old[1]
    return (snew-sold)*(cos(phi-spin1[1])*sin(spin1[2]) + cos(phi-spin2[1])*sin(spin2[2]) + cos(phi-spin3[1])*sin(spin3[2]) + cos(phi-spin4[1])*sin(spin4[2])) + (cnew-cold)*(cos(spin1[2])+cos(spin2[2])+cos(spin3[2])+cos(spin4[2]))   +   d*(cold^2-cnew^2)
end

@inline function EnergyDiffPhi(spin0old::Vector{Float32}, phi::Float32, spin1::Vector{Float32}, spin2::Vector{Float32}, spin3::Vector{Float32}, spin4::Vector{Float32}) #  Like EnergyDiff but optimized for a change in phi only (for XY-phi sweep)
    phiold = spin0old[1]
    return sin(spin0old[2])*((cos(phi-spin1[1])-cos(phiold-spin1[1]))*sin(spin1[2])  +  (cos(phi-spin2[1])-cos(phiold-spin2[1]))*sin(spin2[2])  +  (cos(phi-spin3[1])-cos(phiold-spin3[1]))*sin(spin3[2])  +  (cos(phi-spin4[1])-cos(phiold-spin4[1]))*sin(spin4[2]))
end

function SingleFlip(Lattice::Array{Float32, 3}, i::Int, j::Int, T::Float32, L::Int, σ::Float32, σϕ::Float32, σθ::Float32, pi32::Float32, PBC::Bool, d::Float32, pIsing::Float32, pXY::Float32)      # Try a new configuration (proposal) and either accept or reject it
    newspin, Whichflip = NewSpin(Lattice[:,i,j], σ, σϕ, σθ, pi32, pIsing, pXY)
    if PBC==false
        Δ = 0f0
        if i!=1
            Δ += Correlation(newspin[1], newspin[2], Lattice[1,i-1,j], Lattice[2,i-1,j]) - Correlation(Lattice[1,i,j],Lattice[2,i,j],Lattice[1,i-1,j],Lattice[2,i-1,j]) # This is the right sign as the bounding energy is negative (and is calculate positively)
        end
        if i != L
            Δ += Correlation(newspin[1], newspin[2],Lattice[1,i+1,j],Lattice[2,i+1,j]) - Correlation(Lattice[1,i,j],Lattice[2,i,j],Lattice[1,i+1,j],Lattice[2,i+1,j])
        end
        if j != 1
            Δ += Correlation(newspin[1], newspin[2],Lattice[1,i,j-1],Lattice[2,i,j-1]) - Correlation(Lattice[1,i,j],Lattice[2,i,j],Lattice[1,i,j-1],Lattice[2,i,j-1])
        end
        if j != L
            Δ += Correlation(newspin[1], newspin[2],Lattice[1,i,j+1],Lattice[2,i,j+1]) - Correlation(Lattice[1,i,j],Lattice[2,i,j],Lattice[1,i,j+1],Lattice[2,i,j+1])
        end
        Δ += d * (cos(Lattice[2,i,j])^2 - cos(newspin[2])^2)
    else
        if Whichflip==1
            Δ = EnergyDiff(Lattice[:,i,j], newspin, Lattice[:,mod(i-2,L)+1,j],Lattice[:,i,mod(j-2,L)+1],Lattice[:,mod(i,L)+1,j],Lattice[:,i,mod(j,L)+1],d)      # previously: Δ = Correlation(newspin,Lattice[:,mod(i-2,L)+1,j]) - Correlation(Lattice[:,i,j],Lattice[:,mod(i-2,L)+1,j])   +   Correlation(newspin,Lattice[:,mod(i,L)+1,j]) - Correlation(Lattice[:,i,j],Lattice[:,mod(i,L)+1,j])   +   Correlation(newspin,Lattice[:,i,mod(j-2,L)+1]) - Correlation(Lattice[:,i,j],Lattice[:,i,mod(j-2,L)+1])   +   Correlation(newspin,Lattice[:,i,mod(j,L)+1]) - Correlation(Lattice[:,i,j],Lattice[:,i,mod(j,L)+1])
        elseif Whichflip==3
            Δ = EnergyDiffPhi(Lattice[:,i,j], newspin[1], Lattice[:,mod(i-2,L)+1,j],Lattice[:,i,mod(j-2,L)+1],Lattice[:,mod(i,L)+1,j],Lattice[:,i,mod(j,L)+1])
        else # if Ising or XY-theta flip, as on both cases phi stay the same
            Δ = EnergyDiffTheta(Lattice[:,i,j], newspin[2], Lattice[:,mod(i-2,L)+1,j],Lattice[:,i,mod(j-2,L)+1],Lattice[:,mod(i,L)+1,j],Lattice[:,i,mod(j,L)+1],d)
        end
    end
    expdiff = exp(Δ/T)
    if expdiff > 1 || expdiff > rand(TaskLocalRNG())    # the second condition would be enough to make it works, but this is faster as sometimes it doesn't have to generate a random number
        Lattice[1,i,j] = mod(newspin[1], 2*pi32)
        Lattice[2,i,j] = newspin[2]
        acceptance=true
    else
        acceptance=false
    end
    return Lattice[:,i,j], acceptance, Whichflip
end

function NewSpin(spin::Vector{Float32}, σ::Float32, σϕ::Float32, σθ::Float32, pi32::Float32, pIsing::Union{Float32, Int}, pXY::Union{Float32, Int}) # Choose which kind of proposal (isotropic one, Ising flip, only change phi, only change theta)
    if pIsing == 0 && pXY == 0
        Whichflip = 1
        spin = NewPhiTheta(spin, σ, pi32)
    else
        random = rand(TaskLocalRNG())       # This way to generate less random number, thus faster code
        if pIsing > 0 && random < pIsing
            Whichflip = 2
            spin[2] = pi32 - spin[2]
        elseif random < pXY
            if random > .3*pXY
                spin[1] = NewPhi(spin[1], σϕ, pi32)
                Whichflip = 3
            else
                spin[2] = NewTheta(spin[2], σθ, pi32)
                Whichflip = 4
            end
        else                                # Same thing than the first part of the function, to lose less time (as generate a random number takes way more time than an if / else
            Whichflip = 1
            spin = NewPhiTheta(spin, σ, pi32)
        end
    end
    return spin, Whichflip
end

@inline function NewPhi(phi::Float32, σϕ::Float32, pi32::Float32) # Propose a new phi
    mod(rand(TaskLocalRNG(), Normal(phi, σϕ)), 2*pi32)
end

@inline function NewTheta(theta::Float32, σθ::Float32, pi32::Float32) # Propose a new Theta
    theta = mod(rand(TaskLocalRNG(), Normal(theta, σθ)), 2*pi32)
    if theta > pi32; theta = 2*pi32 - theta
    end
    return theta
end

function NewPhiTheta(spin::Vector{Float32}, σ::Float32, pi32::Float32) # Propose new phi and theta (isotropic sweep)
    rng = TaskLocalRNG()
    alpha  = rand(rng, Float32)*2*pi32
    @fastmath begin
        beta   = min(abs(rand(rng, Normal(0f0, σ))), pi32)
        s1, c1 = sincos(spin[1])
        s2, c2 = sincos(spin[2])
        sa, ca = sincos(alpha)
        sb, cb = sincos(beta)
        X = c2*sb*ca + s2*cb
        spin[1] = atan(s1*X + c1*sb*sa,   c1*X - s1*sb*sa) # alpha is phi'
        spin[2] = acos(clamp(c2*cb - s2*sb*ca, -1f0, 1f0))
    end
    return spin
end

function LatticeSweep(Lattice::Array{Float32, 3}, L::Int, E::Float32, T::Float32, σ::Float32, σϕ::Float32, σθ::Float32, pi32::Float32, PBC::Bool, d::Float32, pIsing::Float32, pXY::Float32)    #MH lattice sweep
    Accept, Tr = zeros(Int, 4), zeros(Int, 4)
    for j in 1:L, k in 1:L
        Lattice[:,j,k], a, Whichflip = SingleFlip(Lattice,j,k,T,L,σ, σϕ, σθ, pi32, PBC, d, pIsing, pXY)
        Accept[Whichflip]+=a
        Tr[Whichflip]+=1
    end
    return Lattice, Accept, Tr
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

function AddNeighbors(i::Int,j::Int, L::Int, InCluster::BitMatrix, p::Float32, ZLattice::Array{Float32, 2}, pi32::Float32)    # For Wolff update: add or not to the cluster the neighbors of a site already in the cluster
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

function Swap(Replicas::Vector{Replica}, i::Int)    # For parallel tempering: try to swap with the next temperature
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

function MHStep(step::Int, rep::Replica, L::Int, d::Float32, PBC::Bool, pi32::Float32)   # ful MH step: i.e. lattice sweep or Wolff algorithm, and (potentially) adjust the three σ for the three gaussian proposal to be closer to a 50% acceptance rate
    Mz2 = mean(cos.(rep.Lattice[2,:,:]).^2)
    if mod(step,1000)==0 && Mz2 > .7f0
        rep.Lattice, rep.E = MultipleIsingFlips(rep.Lattice, L, rep.T, pi32, d, PBC)    # Wolff algorithm
    else
        rep.Lattice, Accept, Tr = LatticeSweep(rep.Lattice, L, rep.E, rep.T, rep.σ, rep.σϕ, rep.σθ, pi32, PBC, d, IsingProba(Mz2), XYProba(Mz2))
        if Tr[1] > 10
            rep.σ = min(max(rep.σ*2*Accept[1]/Tr[1], 1/sqrt(step+1000)), 10)     # if the acceptance is higher than .5, sigma should increase (decreased) to explore more (less) the phasespace to get closer to .5
        end
        if Tr[3] > 0
            a = 1/(Tr[3]+2)
            rep.σϕ = min(max(rep.σϕ*2*min(max(Accept[3]/Tr[3], a), 1-a), 1/sqrt(step + 1000)), 10)     # if the acceptance is higher than .5, sigma should increase (decreased) to explore more (less) the phasespace to get closer to .5
        end         # I could have used clamp(a,b,c) instead of min(max(a,c),b) but it is less efficient (around 50% slower)
        if Tr[4] > 0
            a = 1/(Tr[4]+2)
            rep.σθ = min(max(rep.σθ*2*min(max(Accept[4]/Tr[4], a), 1-a), 1/sqrt(step+10000)), 2)     # if the acceptance is higher than .5, sigma should increase (decreased) to explore more (less) the phasespace to get closer to .5
        end
        rep.Acceptance += Accept
        rep.Try += Tr
    end
    return rep
end

function MH_parallel_tempering(L::Int, N::Int, T::Vector{Float32}, burn::Int, pi32::Float32, AllLattices::Vector{Int}, d::Float32, Save::Bool=false, Skip::Int=10, swap::Int=101)     # Sampler
    nT = length(T)
    # Lattices = Array{Array{Float64,3}}(undef, nT, round(Int,2000/Skip))#-1)       # To save the last lattices
    Replicas = Vector{Replica}(undef, nT)
    for i=1:nT                              # Initialize Replicas
        lat = Initial_lattice(L, pi32)
        Replicas[i] = Replica(lat, T[i], Energy(lat, L, PBC, d), .25f0, .4f0, .25f0, [0,0,0,0], [0,0,0,0])
    end
    Nswap = round(Int, N/swap)                      # number of lattice swap's try 
    Nmeasurement = round(Int, (N - burn)/Skip)      # Number of measurement (as some sweeps are not measured to save some time)
    Energies = zeros(Float32, nT, Nmeasurement)
    Mag  = zeros(Float32, nT, Nmeasurement)
    corr = [zeros(Float32, length(RowCorr(Initial_lattice(L, pi32), L, PBC))) for _ in 1:nT]
    SwapAcceptance = zeros(Int, nT-1)
    
    for i=1:Nswap
        Threads.@threads for n=1:nT
            for j=1:swap
                k = (i-1)*swap+j
                Replicas[n] = MHStep(k, Replicas[n], L, d, PBC, pi32)
                if mod(k, Skip) == 0 && k > burn            # to not measure every lattice sweeps
                    lat = Replicas[n].Lattice
                    m = round(Int,(k-burn)/Skip)
                    Energies[n, m] = Energy(lat, L, PBC, d)
                    Mag[n, m] = sqrt(mean(cos(phi)*sin(theta) for phi in lat[1,:,:], theta in lat[2,:,:])^2 + mean(sin(phi)*sin(theta) for phi in lat[1,:,:], theta in lat[2,:,:])^2 + mean(cos(theta) for theta in lat[2,:,:])^2)
                    corr[n] += RowCorr(lat, L, PBC)
                    # if k > N-2000
                    #     kk = round(Int,(k+2000-N)/Skip)
                    #     Lattices[n, kk] = lat
                    # end
                end
            end
            Replicas[n].E = Energy(Replicas[n].Lattice, L, PBC, d)  #   Only need to calculate the energy for storing (above) and for potentially swapping lattices
        end
        for n=1:nT-1
            SwapAcceptance[n] += Swap(Replicas, n)
        end
    end
    SwapAccept = SwapAcceptance/Nswap
    Energies /= !PBC*L*(L-1) + PBC*L^2
    corr /= Nmeasurement
    if Save==true
        for n=1:nT
            accept = Replicas[n].Acceptance ./ Replicas[n].Try
            @save "Data/$(AllLattices)_$N/$(L)_$(T[n])_$d.jld2" Energies=Energies[n,:] Mag=Mag[n,:] corr=corr[n] accept# Lattices=Lattices[n,:]# Energies2
        end
        @save "Data/$(AllLattices)_$N/swap_$(L)_$d.jld2" SwapAccept
    end
    # for i in eachindex(Lattices)
    #     if !isassigned(Lattices,i); println(i);end
    # end
    for n=1:nT; println("$L \t $d \t $(T[n]) \t $(Replicas[n].Acceptance ./ Replicas[n].Try) \t and Magz2 : $(mean(cos.(Replicas[n].Lattice[2,:,:]).^2))"); end
    return mean(Energies[1,:]), #=mean(Energies2[1,:]),=# mean(Mag[1,:])
end

function RowCorr(Lattice::Array{Float32, 3}, L::Int, PBC::Bool)   # Return the average row Correlation of the matrix (in a Correlation vector : first neigbour, second neigbour...)
    if PBC==false
        N = L-1
        corr = zeros(Float32, N)
        @inbounds for i=1:N
            a=0f0
            for j=1:L
                for k=1:(L-i)
                    a+=Correlation(Lattice[1,j,k], Lattice[2,j,k], Lattice[1,j,k+i], Lattice[2,j,k+i])
                end
            end
            corr[i] = a/((L-i)*L)
        end
    else
        N=round(Int, L/2-1)
        corr = zeros(Float32, N)
        @inbounds for i=1:N
            a=0f0
            for j=1:L
                for k=1:L
                    k2 = (k + i <= L) ? k + i : k + i - L # like modulo, but more efficient (only if x < 2y for mod(x,y))
                    a+=Correlation(Lattice[1,j,k], Lattice[2,j,k], Lattice[1,j,k2], Lattice[2,j,k2]) #cos(Lattice[1,j,k]-Lattice[1,j,mod(k+i-1,L)+1]) * cos(Lattice[2,j,k]-Lattice[2,j,mod(k+i-1,L)+1])
                end
            end
            corr[i] = a/L^2
        end
    end
    return corr
end

@inline function Temperatures(Tmin::Float64, Tmax::Float64, N::Int)
    return Float32.(round.(Tmin .* (Tmax/Tmin) .^ ((0:N-1)/(N-1)); digits=3))
end

@inline @fastmath function IsingProba(Mz2::Float32)   # Probability of Ising-flip, just a choice, there is no perfect one
    return (1.9f0*Mz2 .-.95f0).^2 .*(sign.(Mz2 .- .5f0).+1)/2
end

@inline @fastmath function XYProba(Mz2::Float32)      # Probability of XY-flip (either phi or theta, independently): useful to converge faster for XY config as σϕ and σθ are tuned independently, to really have θ very close to π/2
    return (.9f0 .- 9f0*Mz2).^2 .*(sign.(.1f0 .-Mz2).+1)/2
end # just a choice, there is no perfect one



#   @inline     directly put the code of the function where the function is called, instead of really calling it    for single line function, no need "return" inside it
#   @fastmath   allows to reorder operation to save a bit of time (sometimes less precise because of floating point handling)
#   @inbounds   don't check if wanted component exist example     x=vector[n]     Doesn't check that the n-th component even exists, just trust the process