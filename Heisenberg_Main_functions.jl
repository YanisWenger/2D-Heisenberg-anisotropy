mutable struct Replica
    Lattice::Array{Float32,3}
    T::Float32                  # Constant
    E::Float32                  # Instantaneous
    σ::Float32                  # Instantaneous
    σϕ::Float32                 # Instantaneous
    σθ::Float32                 # Instantaneous
    Acceptance::Vector{Int64}   # History
    Try::Vector{Int64}          # History
    Cluster::Int64              # History (number of cluster update)
end

function Energy(Lattice::Array{Float32, 3}, L::Int, PBC::Bool, D::Float32)    # Calculate the Energy of a configuration
    E=0f0
    n = PBC ? L : L-1
    @inbounds for i=1:n
        i2 = i==L ? 1 : i+1
        @simd for j=1:n
            j2 = j==L ? 1 : j+1
            Lat_ij = @SVector [Lattice[1, i, j], Lattice[2, i, j]]
            neighbor_i = @SVector [Lattice[1, i, j2], Lattice[2, i, j2]]
            neighbor_j = @SVector [Lattice[1, i2, j], Lattice[2, i2, j]]
            E += -Correlation(Lat_ij, neighbor_i) - Correlation(Lat_ij, neighbor_j) + D * cos(Lat_ij[2])^2
        end
    end
    return E
end #   @inbounds almost divide the execution time by 2

@inline @fastmath @inbounds function Correlation(spin1::SVector{2, Float32}, spin2::SVector{2, Float32})
    sθ1, cθ1 = sincos(spin1[2])
    sθ2, cθ2 = sincos(spin2[2])
    return sθ1*sθ2*cos(spin1[1]-spin2[1]) + cθ1*cθ2
end

function Initial_lattice(L::Int, rng::AbstractRNG)
    lattice = Array{Float32}(undef, 2, L, L)
    lattice[1, :, :] .= 2*Float32(π)*rand(rng, Float32, L, L)
    lattice[2, :, :] .= acos.(1 .- 2*rand(rng, Float32, L, L))
    return lattice
end

@inline @fastmath @inbounds function EnergyDiff(spin0old::SVector{2, Float32}, spin0new::SVector{2, Float32}, spin1::SVector{2, Float32}, spin2::SVector{2, Float32}, spin3::SVector{2, Float32}, spin4::SVector{2, Float32}, D::Float32) # calculate the bounding energy difference between a spin and its neigbors and another one with the same neighbors
    sold, snew = sin(spin0old[2]), sin(spin0new[2])
    cold, cnew = cos(spin0old[2]), cos(spin0new[2])
    return (snew*cos(spin0new[1]-spin1[1]) - sold*cos(spin0old[1]-spin1[1]))*sin(spin1[2]) + (cnew-cold)*cos(spin1[2])  +  (snew*cos(spin0new[1]-spin2[1]) - sold*cos(spin0old[1]-spin2[1]))*sin(spin2[2]) + (cnew-cold)*cos(spin2[2])  +  (snew*cos(spin0new[1]-spin3[1]) - sold*cos(spin0old[1]-spin3[1]))*sin(spin3[2]) + (cnew-cold)*cos(spin3[2])  +  (snew*cos(spin0new[1]-spin4[1]) - sold*cos(spin0old[1]-spin4[1]))*sin(spin4[2]) + (cnew-cold)*cos(spin4[2])   +   D*(cold^2-cnew^2) # old - new for the last term, as H = - Σ Si Sj + D Σ Sz^2 (sign difference between these two terms)
end

@inline @fastmath @inbounds function EnergyDiffTheta(spin0old::SVector{2, Float32}, theta::Float32, spin1::SVector{2, Float32}, spin2::SVector{2, Float32}, spin3::SVector{2, Float32}, spin4::SVector{2, Float32}, D::Float32) # Like EnergyDiff but optimized for a change in theta only (for XY-theta sweep or Ising flip)
    sold, snew = sin(spin0old[2]), sin(theta)
    cold, cnew = cos(spin0old[2]), cos(theta)
    phi = spin0old[1]
    return (snew-sold)*(cos(phi-spin1[1])*sin(spin1[2]) + cos(phi-spin2[1])*sin(spin2[2]) + cos(phi-spin3[1])*sin(spin3[2]) + cos(phi-spin4[1])*sin(spin4[2])) + (cnew-cold)*(cos(spin1[2])+cos(spin2[2])+cos(spin3[2])+cos(spin4[2]))   +   D*(cold^2-cnew^2)
end

@inline @fastmath @inbounds function EnergyDiffPhi(spin0old::SVector{2, Float32}, phi::Float32, spin1::SVector{2, Float32}, spin2::SVector{2, Float32}, spin3::SVector{2, Float32}, spin4::SVector{2, Float32}) #  Like EnergyDiff but optimized for a change in phi only (for XY-phi sweep)
    phiold = spin0old[1]
    return sin(spin0old[2])*((cos(phi-spin1[1])-cos(phiold-spin1[1]))*sin(spin1[2])  +  (cos(phi-spin2[1])-cos(phiold-spin2[1]))*sin(spin2[2])  +  (cos(phi-spin3[1])-cos(phiold-spin3[1]))*sin(spin3[2])  +  (cos(phi-spin4[1])-cos(phiold-spin4[1]))*sin(spin4[2]))
end

@inbounds function SingleFlip(Lattice::Array{Float32, 3}, i::Int, j::Int, T::Float32, L::Int, σ::Float32, σϕ::Float32, σθ::Float32, PBC::Bool, D::Float32, pIsing::Float32, pXY::Float32, rng::AbstractRNG)      # Try a new configuration (proposal) and either accept or reject it
    newspin, Whichflip = NewSpin(Lattice[:,i,j], σ, σϕ, σθ, pIsing, pXY, rng)
    if PBC==false
        Lat_ij = @SVector [Lattice[1,i,j], Lattice[2,i,j]]
        Δ = 0f0
        if i!=1
            spin2 = @SVector [Lattice[1,i-1,j], Lattice[2,i-1,j]]
            Δ += Correlation(newspin, spin2) - Correlation(Lat_ij,spin2) # This is the right sign as the bounding energy is negative (and is calculate positively)
        end
        if i != L
            spin2 = @SVector [Lattice[1,i+1,j], Lattice[2,i+1,j]]
            Δ += Correlation(newspin, spin2) - Correlation(Lat_ij,spin2)
        end
        if j != 1
            spin2 = @SVector [Lattice[1,i,j-1], Lattice[2,i,j-1]]
            Δ += Correlation(newspin, spin2) - Correlation(Lat_ij,spin2)
        end
        if j != L
            spin2 = @SVector [Lattice[1,i,j+1], Lattice[2,i,j+1]]
            Δ += Correlation(newspin, spin2) - Correlation(Lat_ij,spin2)
        end
        Δ += D * (cos(Lattice[2,i,j])^2 - cos(newspin[2])^2)
    else
        Lat_ij = @SVector [Lattice[1,i,j], Lattice[2,i,j]]
        neighbor1 = @SVector [Lattice[1,i==1 ? L : i-1,j], Lattice[2,i==1 ? L : i-1,j]]
        neighbor2 = @SVector [Lattice[1,i,j==1 ? L : j-1], Lattice[2,i,j==1 ? L : j-1]]
        neighbor3 = @SVector [Lattice[1,i==L ? 1 : i+1,j], Lattice[2,i==L ? 1 : i+1,j]]
        neighbor4 = @SVector [Lattice[1,i,j==L ? 1 : j+1], Lattice[2,i,j==L ? 1 : j+1]]
        if Whichflip==1
            Δ = EnergyDiff(Lat_ij, newspin, neighbor1, neighbor2, neighbor3, neighbor4,D)      # previously: Δ = Correlation(newspin,Lattice[:,mod(i-2,L)+1,j]) - Correlation(Lattice[:,i,j],Lattice[:,mod(i-2,L)+1,j])   +   Correlation(newspin,Lattice[:,mod(i,L)+1,j]) - Correlation(Lattice[:,i,j],Lattice[:,mod(i,L)+1,j])   +   Correlation(newspin,Lattice[:,i,mod(j-2,L)+1]) - Correlation(Lattice[:,i,j],Lattice[:,i,mod(j-2,L)+1])   +   Correlation(newspin,Lattice[:,i,mod(j,L)+1]) - Correlation(Lattice[:,i,j],Lattice[:,i,mod(j,L)+1])
        elseif Whichflip==3
            Δ = EnergyDiffPhi(Lat_ij, newspin[1], neighbor1, neighbor2, neighbor3, neighbor4)
        else # if Ising or XY-theta flip, as on both cases phi stay the same
            Δ = EnergyDiffTheta(Lat_ij, newspin[2], neighbor1, neighbor2, neighbor3, neighbor4,D)
        end
    end
    expdiff = exp(Δ/T)
    if expdiff > 1 || expdiff > rand(rng)    # the second condition would be enough to make it works, but this is faster as sometimes it doesn't have to generate a random number
        Lattice[1,i,j] = mod(newspin[1], 2*Float32(π))
        Lattice[2,i,j] = newspin[2]
        acceptance=true
    else
        acceptance=false
    end
    return Lattice[:,i,j], acceptance, Whichflip
end

@inbounds function NewSpin(spin::Vector{Float32}, σ::Float32, σϕ::Float32, σθ::Float32, pIsing::Union{Float32, Int}, pXY::Union{Float32, Int}, rng::AbstractRNG) # Choose which kind of proposal (isotropic one, Ising flip, only change phi, only change theta)
    if pIsing == 0 && pXY == 0
        Whichflip = 1
        newspin = NewPhiTheta(spin, σ, rng)
    else
        random = rand(rng)       # This way to generate less random number, thus faster code
        if pIsing > 0 && random < pIsing
            Whichflip = 2
            newspin = @SVector [spin[1], Float32(π) - spin[2]]
        elseif random < pXY
            if random > .3*pXY
                newspin = @SVector [NewPhi(spin[1], σϕ, rng), spin[2]]
                Whichflip = 3
            else
                newspin = @SVector [spin[1], NewTheta(spin[2], σθ, rng)]
                Whichflip = 4
            end
        else                    # Same thing than the first part of the function, to lose less time (as generate a random number takes way more time than an if / else
            Whichflip = 1
            newspin = NewPhiTheta(spin, σ, rng)
        end
    end
    return newspin, Whichflip
end

@inline function NewPhi(phi::Float32, σϕ::Float32, rng::AbstractRNG) # Propose a new phi
    mod(rand(rng, Normal(phi, σϕ)), 2*Float32(π))
end

@inline function NewTheta(theta::Float32, σθ::Float32, rng::AbstractRNG) # Propose a new Theta
    pi32 = Float32(π)
    theta = mod(rand(rng, Normal(theta, σθ)), 2*pi32)
    if theta > pi32; theta = 2*pi32 - theta
    end
    return theta
end

function NewPhiTheta(spin::Vector{Float32}, σ::Float32, rng::AbstractRNG) # Propose new phi and theta (isotropic sweep)
    pi32 = Float32(π)
    alpha  = rand(rng, Float32)*2*pi32
    @fastmath begin
        beta   = min(abs(rand(rng, Normal(0f0, σ))), pi32)
        s1, c1 = sincos(spin[1])
        s2, c2 = sincos(spin[2])
        sa, ca = sincos(alpha)
        sb, cb = sincos(beta)
        X = c2*sb*ca + s2*cb
        spin = @SVector [atan(s1*X + c1*sb*sa,   c1*X - s1*sb*sa), acos(clamp(c2*cb - s2*sb*ca, -1f0, 1f0))] # alpha is phi'
    end
    return spin
end

function LatticeSweep(Lattice::Array{Float32, 3}, L::Int, T::Float32, σ::Float32, σϕ::Float32, σθ::Float32, PBC::Bool, d::Float32, pIsing::Float32, pXY::Float32, rng::AbstractRNG)    #MH lattice sweep
    Accept, Tr = zeros(Int, 4), zeros(Int, 4)
    @inbounds for k in 1:L, j in 1:L # faster as Julia stores arrays in column-major order
        Lattice[:,j,k], a, Whichflip = SingleFlip(Lattice,j,k,T,L,σ, σϕ, σθ, PBC, d, pIsing, pXY, rng)
        Accept[Whichflip]+=a
        Tr[Whichflip]+=1
    end
    return Lattice, Accept, Tr
end

function IsingClusterUpdate(ZLattice::AbstractArray{Float32, 2}, L::Int, T::Float32, rng::AbstractRNG)      # Wolff update
    pi32 = Float32(π)
    β2 = 2/T
    k = rand(rng, 1:L)
    l = rand(rng, 1:L)
    InCluster = falses(L,L)
    InCluster[k,l] = true
    queue = Vector{Tuple{Int, Int}}(undef, L*L) # which site to check neighbor of
    Head = 1    # site that is being checked
    Tail = 1    # site to check after
    queue[Tail] = (k,l)
    Tail += 1
    @inbounds while Head < Tail
        i, j = queue[Head]
        Head += 1
        neighbors = ((i==1 ? L : i-1,j), (i,j==1 ? L : j-1), (i==L ? 1 : i+1,j), (i,j==L ? 1 : j+1))    # The four neighbors (could be optimized using a?b:c instead of mod)
        cos_site = cos(ZLattice[i,j])
        for n=1:4               # To check the four neighbors
            ki,kj = neighbors[n]
            if !InCluster[ki,kj] && rand(rng) < 1-exp(-β2*cos_site*cos(ZLattice[ki,kj]))
                InCluster[ki,kj] = true
                queue[Tail] = (ki, kj)  # adding site to the queue
                Tail += 1
            end
        end
    end
    @inbounds for j in 1:L, i in 1:L # faster as Julia stores arrays in column-major order
        if InCluster[i,j]
            ZLattice[i,j] = pi32 - ZLattice[i,j]
        end
    end
    return ZLattice
end

# function ClusterUpdateIsotropic(Lattice::Array{Float32, 3}, L::Int, T::Float32, rng::AbstractRNG)      # Wolff update
#     β2 = 2/T
#     ϕ = 2*Float32(π)*rand(rng, Float32)     # random Axis
#     θ = acos.(1 .- 2*rand(rng, Float32))
#     k = rand(rng, 1:L)
#     l = rand(rng, 1:L)
#     # ϕ = Lattice[1,k,l]                    # Axis aligned with the first spin of the cluster (chosen randomly)
#     # θ = Lattice[2,k,l]
#     Axis = @SVector [ϕ, θ]    # random direction to flip the cluster
#     InCluster = falses(L,L)
#     InCluster[k,l] = true
#     queue = Vector{Tuple{Int, Int}}(undef, L*L) # which site to check neighbor of
#     Head = 1    # site that is being checked
#     Tail = 1    # site to check after
#     queue[Tail] = (k,l)
#     Tail += 1
#     while Head < Tail
#         i, j = queue[Head]
#         Head += 1
#         neighbors = ((i==1 ? L : i-1,j), (i,j==1 ? L : j-1), (i==L ? 1 : i+1,j), (i,j==L ? 1 : j+1))    # The four neighbors (could be optimized using a?b:c instead of mod)
#         dot0 = Correlation(Axis, @SVector [Lattice[1,i,j], Lattice[2,i,j]])
#         for n=1:4               # To check the four neighbors
#             ki,kj = neighbors[n]
#             if !InCluster[ki,kj] && rand(rng) < 1-exp(-β2*dot0*Correlation(Axis, @SVector [Lattice[1,ki,kj],Lattice[2,ki,kj]]))
#                 InCluster[ki,kj] = true
#                 queue[Tail] = (ki, kj)  # adding site to the queue
#                 Tail += 1
#             end
#         end
#     end
#     twoθ = 2*θ
#     twoϕ = 2*ϕ
#     @inbounds for i in 1:L, j in 1:L
#         if InCluster[i,j]
#             Lattice[2,i,j], turnPhi = FlipTheta(twoθ, Lattice[2,i,j])
#             Lattice[1,i,j] = FlipPhi(twoϕ, Lattice[1,i,j], turnPhi)
#         end
#     end
#     return Lattice
# end

function FlipTheta(twoθ::Float32, ThetaSpin::Float32)
    pi32 = Float32(π)
    x = twoθ - ThetaSpin
    if x > 0 && x < pi32
        ThetaSpin = pi32 - x
        turnPhi = true
    else
        if x <= 0
            ThetaSpin = pi32 + x
        else    # x can only be between -pi and 2pi, thus, this is equivalent to if pi < x < 2pi
            ThetaSpin = x - pi32
        end
        turnPhi = false
    end
    return ThetaSpin, turnPhi
end

function FlipPhi(twoϕ::Float32, PhiSpin::Float32, turnPhi::Bool)
    pi32 = Float32(π)
    x = twoϕ - PhiSpin + pi32*(1+turnPhi) # between -pi and 6pi, thus more efficient to do if elseif than modulo
    if x < 0
        PhiSpin = x + 2*pi32
    elseif x < 2*pi32
        PhiSpin = x
    elseif x < 4*pi32
        PhiSpin = x - 2*pi32
    else
        PhiSpin = x - 4*pi32
    end
    return PhiSpin
end

function Swap(Replicas::Vector{Replica}, i::Int, rng_swap::AbstractRNG, βdiff_i::Float32)    # For parallel tempering: try to swap with the next temperature
    j=i+1
    if exp(βdiff_i*(Replicas[j].E - Replicas[i].E)) > rand(rng_swap)
        Replicas[i].Lattice, Replicas[j].Lattice = Replicas[j].Lattice, Replicas[i].Lattice
        Replicas[i].E,       Replicas[j].E       = Replicas[j].E,       Replicas[i].E
        return true
    else
        return false
    end
end

function MHStep(step::Int, rep::Replica, L::Int, D::Float32, PBC::Bool, rng::AbstractRNG)   # ful MH step: i.e. lattice sweep or Wolff algorithm, and (potentially) adjust the three σ for the three gaussian proposal to be closer to a 50% acceptance rate
    Mz2 = Mag_z2(rep.Lattice, L)
    IsingProba_Mz2 = IsingProba(Mz2)
    if mod(step,100)==0 && IsingProba_Mz2 > 0 && IsingProba_Mz2 > rand(rng) # the second condition is mathematically useless because of the third, but it doesn't generate random number if proba =0
        rep.Lattice[2,:,:] = IsingClusterUpdate(view(rep.Lattice, 2, : , : ), L, rep.T, rng)    # Wolff algorithm
        rep.Cluster += 1
    else
        rep.Lattice, Accept, Tr = LatticeSweep(rep.Lattice, L, rep.T, rep.σ, rep.σϕ, rep.σθ, PBC, D, IsingProba_Mz2, XYProba(Mz2), rng)
        @fastmath begin
            if Tr[1] > 10
                rep.σ = min(max(rep.σ*2.86*Accept[1]/Tr[1], 1/sqrt(step+1000)), 10)     # if the acceptance is higher (lower) than .35, sigma should increase (decreased) to explore more (less) the phase space to get closer to .35
            end
            if Tr[3] > 0
                a = 1/(Tr[3]+2)
                rep.σϕ = min(max(rep.σϕ*2.27*min(max(Accept[3]/Tr[3], a), 1-a), 1/sqrt(step + 1000)), 10)     # if the acceptance is higher (lower) than .44, sigma should increase (decreased) to explore more (less) the phase space to get closer to .44
            end         # I could have used clamp(a,b,c) instead of min(max(a,c),b) but it is less efficient (around 50% slower)
            if Tr[4] > 0
                a = 1/(Tr[4]+2)
                rep.σθ = min(max(rep.σθ*2.27*min(max(Accept[4]/Tr[4], a), 1-a), 1/sqrt(step+10000)), 2)     # if the acceptance is higher (lower) than .44, sigma should increase (decreased) to explore more (less) the phase space to get closer to .44
            end
            rep.Acceptance += Accept
            rep.Try += Tr
        end
    end
    return rep
end

function MH_parallel_tempering(L::Int, N::Int, T::Vector{Float32}, burn::Int, AllLattices::Vector{Int}, D::Float32, Save::Bool=false, SaveLattices::Bool=false, Skip::Int=10, swap::Int=80)     # Sampler
    nT = length(T)
    Lattices = Array{Array{Float64,3}, 2}(undef, nT, div(2000, Skip)-!(N%swap==0)*Skip)       # To save the last lattices, 
    Replicas = Vector{Replica}(undef, nT)
    βdiff = diff(1 ./T)    # inverse temperature difference, useful to swap lattices
    seed = rand(1:10000)
    rngs = [Xoshiro(seed+i) for i in 1:nT]
    Threads.@threads for i=1:nT                              # Initialize Replicas
        lat = Initial_lattice(L, rngs[i])
        Replicas[i] = Replica(lat, T[i], Energy(lat, L, PBC, D), .25f0, .4f0, .25f0, [0,0,0,0], [0,0,0,0], 0)
    end
    Nswap = div(N, swap)                            # number of lattice swap's try 
    println("$N, $Nswap, $swap, $(Nswap*swap)")
    Nmeasurement = div(N - burn, Skip)      # Number of measurement (as some sweeps are not measured to save some time)
    Energies = Array{Float32}(undef, nT, Nmeasurement)
    Mag  = Array{Float32}(undef, nT, Nmeasurement, 3)
    corr = [zeros(Float32, length(RowCorr(Replicas[1].Lattice, L, PBC))) for _ in 1:nT]
    SwapAcceptance = zeros(Int, nT-1)
    
    for i=1:Nswap
        Threads.@threads for n=1:nT
            for j=1:swap
                k = (i-1)*swap+j
                Replicas[n] = MHStep(k, Replicas[n], L, D, PBC, rngs[n])
                if k > burn && mod(k, Skip) == 0            # to not measure every lattice sweeps
                    lat = Replicas[n].Lattice
                    m = div(k-burn, Skip)
                    Energies[n, m] = Energy(lat, L, PBC, D)
                    Mag[n, m, :] = mag(lat, L)
                    corr[n] += RowCorr(lat, L, PBC)
                    if SaveLattices && k > N-2000
                        kk = div(k+2000-N, Skip)
                        Lattices[n, kk] = lat
                    end
                end
            end
            Replicas[n].E = Energy(Replicas[n].Lattice, L, PBC, D)  #   Only need to calculate the energy for storing (above) and for potentially swapping lattices
        end
        for n=1:nT-1
            SwapAcceptance[n] += Swap(Replicas, n, rngs[1], βdiff[n])
        end
    end
    SwapAccept = SwapAcceptance/Nswap
    Energies /= !PBC*L*(L-1) + PBC*L^2
    corr /= Nmeasurement
    if Save==true
        for n=1:nT
            accept = Replicas[n].Acceptance ./ Replicas[n].Try
            if SaveLattices
                @save "Data/$(AllLattices)_$N/$(L)_$(T[n])_$D.jld2" Energies=Energies[n,:] Mx=Mag[n,:,1] My=Mag[n,:,2] Mz=Mag[n,:,3] corr=corr[n] accept Lattices=Lattices[n,:]
            else
                @save "Data/$(AllLattices)_$N/$(L)_$(T[n])_$D.jld2" Energies=Energies[n,:] Mx=Mag[n,:,1] My=Mag[n,:,2] Mz=Mag[n,:,3] corr=corr[n] accept
            end
        end
        @save "Data/$(AllLattices)_$N/swap_$(L)_$D.jld2" SwapAccept
    end
    for n=1:nT; println("$L \t $D \t $(T[n]) \t $(round.(Replicas[n].Acceptance ./ Replicas[n].Try; digits=3)) \t and Magz2 : $(round(mean(cos.(Replicas[n].Lattice[2,:,:]).^2); digits=3)) \t $(Replicas[n].Cluster)"); end
    return Energies[1,end]
end

function RowCorr(Lattice::Array{Float32, 3}, L::Int, PBC::Bool)   # Return the average row Correlation of the matrix (in a Correlation vector : first neigbour, second neigbour...)
    if PBC==false
        N = L-1
        corr = zeros(Float32, N)
        @inbounds @fastmath for i=1:N
            a=0f0
            for j=1:L
                for k=1:(L-i)
                    spin1 = @SVector [Lattice[1,j,k], Lattice[2,j,k]]
                    spin2 = @SVector [Lattice[1,j,k+i], Lattice[2,j,k+i]]
                    a+=Correlation(spin1, spin2)
                end
            end
            corr[i] = a/((L-i)*L)
        end
    else
        N=round(Int, L/2-1)
        corr = zeros(Float32, N)
        @inbounds @fastmath for i=1:N
            a=0f0
            for j=1:L
                for k=1:L
                    k2 = (k + i <= L) ? k + i : k + i - L # like modulo, but more efficient (only if x < 2y for mod(x,y))
                    spin1 = @SVector [Lattice[1,j,k], Lattice[2,j,k]]
                    spin2 = @SVector [Lattice[1,j,k2], Lattice[2,j,k2]]
                    a += Correlation(spin1, spin2) #cos(Lattice[1,j,k]-Lattice[1,j,mod(k+i-1,L)+1]) * cos(Lattice[2,j,k]-Lattice[2,j,mod(k+i-1,L)+1])
                end
            end
            corr[i] = a/L^2
        end
    end
    return corr
end

function mag(Lattice::Array{Float32,3}, L::Int64)
    X = 0f0
    Y = 0f0
    Z = 0f0
    @inbounds @fastmath for k in 1:L, j in 1:L # faster as Julia stores arrays in column-major order
        stheta, ctheta = sincos(Lattice[2,j,k])
        sphi, cphi = sincos(Lattice[1,j,k])
        x = cphi*stheta
        y = sphi*stheta
        z = ctheta
        X +=x
        Y +=y
        Z +=z
    end
    return [X, Y, Z]/L^2
end

function Mag_z2(Lattice::Array{Float32,3}, L::Int64)
    magz2 = 0f0
    @inbounds @fastmath for j=1:L, k=1:L
        magz2 += cos(Lattice[2,j,k])^2
    end
    return magz2/L^2
end

@inline function Temperatures(Tmin::Real, Tmax::Real, N::Int64)
    return Float32.(round.(Tmin .* (Tmax/Tmin) .^ ((0:N-1)/(N-1)); digits=3))
end

@inline @fastmath function IsingProba(Mz2::Float32)   # Probability of Ising-flip, just a choice, there is no perfect one
    return (1.9f0*Mz2 -.95f0)^2 *(sign(Mz2 - .5f0)+1)/2
end

@inline @fastmath function XYProba(Mz2::Float32)      # Probability of XY-flip (either phi or theta, independently): useful to converge faster for XY config as σϕ and σθ are tuned independently, to really have θ very close to π/2
    return (.9f0 - 9f0*Mz2)^2 *(sign(.1f0 -Mz2)+1)/2
end # just a choice, there is no perfect one



#   @inline     directly put the code of the function where the function is called, instead of really calling it    for single line function, no need "return" inside it
#   @fastmath   allows to reorder operation to save a bit of time (sometimes less precise because of floating point handling)
#   @inbounds   don't check if wanted component exist example     x=vector[n]     Doesn't check that the n-th component even exists, just trust the process
# Lattices=Array{Int64, 2}(undef, 1, 2)
# !isassigned(Lattices,i)