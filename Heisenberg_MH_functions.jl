function Energy(Lattice::Array{SVector{2,Float32},2}, L::Int, D::Float32)    # Calculate the Energy of a configuration
    E=0f0
    @inbounds for i=1:L
        i2 = i==L ? 1 : i+1
        @simd for j=1:L
            Lat_ij = Lattice[i, j]
            E += -Correlation(Lat_ij, Lattice[i, j==L ? 1 : j+1]) - Correlation(Lat_ij, Lattice[i2, j]) + D * cos(Lat_ij[2])^2
        end
    end
    return E
end #   @inbounds almost divide the execution time by 2

@inline @fastmath @inbounds function Correlation(spin1::SVector{2, Float32}, spin2::SVector{2, Float32}) # Calculate the correlation between to spin i.e. dot product
    sθ1, cθ1 = sincos(spin1[2])
    sθ2, cθ2 = sincos(spin2[2])
    return sθ1*sθ2*cos(spin1[1]-spin2[1]) + cθ1*cθ2
end

@inline function Initial_lattice(L::Int, rng::AbstractRNG)  # generate a lattice with random spin orientation
    return SVector{2,Float32}.(2f0 * π .* rand(rng, Float32, L, L),    acos.(1f0 .- 2f0 .* rand(rng, Float32, L, L)))
end

@inline @fastmath @inbounds function EnergyDiff(spin0old::SVector{2, Float32}, spin0new::SVector{2, Float32}, spin1::SVector{2, Float32}, spin2::SVector{2, Float32}, spin3::SVector{2, Float32}, spin4::SVector{2, Float32}, D::Float32) # calculate the bounding energy difference between a spin and its neighbors and the energy difference between a potentially new one with the same neighbors
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

@inbounds function SingleFlip(Lattice::Array{SVector{2,Float32},2}, i::Int, j::Int, T::Float32, L::Int, σ::Float32, σϕ::Float32, σθ::Float32, D::Float32, pIsing::Float32, pXY::Float32, rng::AbstractRNG)   # Try a new spin configuration (proposal) and either accept or reject it
    newspin, Whichflip = NewSpin(Lattice[i,j], σ, σϕ, σθ, pIsing, pXY, rng)
    if Whichflip==1
        Δ = EnergyDiff(Lattice[i,j], newspin, Lattice[i==1 ? L : i-1,j], Lattice[i,j==1 ? L : j-1], Lattice[i==L ? 1 : i+1,j], Lattice[i,j==L ? 1 : j+1],D)
    elseif Whichflip==3
        Δ = EnergyDiffPhi(Lattice[i,j], newspin[1], Lattice[i==1 ? L : i-1,j], Lattice[i,j==1 ? L : j-1], Lattice[i==L ? 1 : i+1,j], Lattice[i,j==L ? 1 : j+1])
    else # if Ising or XY-theta flip, as on both cases phi stay the same
        Δ = EnergyDiffTheta(Lattice[i,j], newspin[2], Lattice[i==1 ? L : i-1,j], Lattice[i,j==1 ? L : j-1], Lattice[i==L ? 1 : i+1,j], Lattice[i,j==L ? 1 : j+1],D)
    end
    expdiff = exp(Δ/T)
    if expdiff > 1 || expdiff > rand(rng)    # the second condition would be enough to make it works, but this is faster as sometimes it doesn't have to generate a random number
        Lattice[i,j] = newspin
        acceptance=true
    else
        acceptance=false
    end
    return Lattice[i,j], acceptance, Whichflip
end

@inbounds function NewSpin(spin::SVector{2, Float32}, σ::Float32, σϕ::Float32, σθ::Float32, pIsing::Union{Float32, Int}, pXY::Union{Float32, Int}, rng::AbstractRNG) # Choose which kind of proposal (isotropic one, Ising flip, only change phi, only change theta)
    if pIsing == 0 && pXY == 0
        Whichflip = 1
        newspin = NewPhiTheta(spin, σ, rng)
    else
        random = rand(rng)       # This way to generate less random number, thus faster code
        if random < pIsing
            Whichflip = 2
            newspin = @SVector [spin[1], Float32(π) - spin[2]]
        elseif random < pXY
            if random < .7*pXY
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

function NewPhiTheta(spin::SVector{2, Float32}, σ::Float32, rng::AbstractRNG) # Propose new phi and theta (isotropic sweep)
    @fastmath begin
        pi32 = Float32(π)
        alpha  = rand(rng, Float32)*2*pi32
        beta   = min(abs(rand(rng, Normal(0f0, σ))), pi32)
        s1, c1 = sincos(spin[1])
        s2, c2 = sincos(spin[2])
        sa, ca = sincos(alpha)
        sb, cb = sincos(beta)
        X = c2*sb*ca + s2*cb
        spin = @SVector [atan(s1*X + c1*sb*sa,   c1*X - s1*sb*sa), acos(min(max(c2*cb - s2*sb*ca, -1f0), 1f0))] # need to clamp it (max and min are faster than clamp), as there are sometimes bad approximations
    end
    return spin
end

function LatticeSweep(Lattice::Array{SVector{2,Float32},2}, L::Int, T::Float32, σ::Float32, σϕ::Float32, σθ::Float32, d::Float32, pIsing::Float32, pXY::Float32, rng::AbstractRNG)    # MH lattice sweep
    Accept, Tr = zeros(Int, 4), zeros(Int, 4)
    @inbounds for k in 1:L, j in 1:L # faster as Julia stores arrays in column-major order
        Lattice[j,k], a, Whichflip = SingleFlip(Lattice,j,k,T,L,σ, σϕ, σθ, d, pIsing, pXY, rng)
        Accept[Whichflip]+=a
        Tr[Whichflip]+=1
    end
    return Lattice, Accept, Tr
end

function IsingClusterUpdate(ZLattice::AbstractArray{Float32, 2}, L::Int, T::Float32, rng::AbstractRNG)      # z-axis Wolff update
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
    @inbounds for j in 1:L, i in 1:L # faster j, i than i, j as Julia stores arrays in column-major order
        if InCluster[i,j]
            ZLattice[i,j] = pi32 - ZLattice[i,j]
        end
    end
    return ZLattice
end

function MHStep(step::Int, lat::Matrix{SVector{2, Float32}}, T::Float32, L::Int, D::Float32, σ::Float32, σϕ::Float32, σθ::Float32, Acceptance::Vector{Int64}, Try::Vector{Int64}, Cluster_update::Int64, rng::AbstractRNG)   # ful MH step: i.e. lattice sweep or Wolff algorithm, and (potentially) adjust the three σ for the three gaussian proposal to be closer to the optimal acceptance rate
    Mz2 = Mag_z2(getindex.(lat, 2), L)  # quantifies how much spins are aligned along z
    IsingProba_Mz2 = IsingProba(Mz2)
    if mod(step,100)==0 && IsingProba_Mz2 > 0 && IsingProba_Mz2 > rand(rng) # the second condition is mathematically useless because of the third, but it is faster like so as it doesn't generate random number if Probability is 0
        lat = SVector{2,Float32}.(getindex.(lat, 1),    IsingClusterUpdate(getindex.(lat, 2), L, T, rng))    # Wolff algorithm
        Cluster_update += 1
    else
        lat, Accept, Tr = LatticeSweep(lat, L, T, σ, σϕ, σθ, D, IsingProba_Mz2, XYProba(Mz2), rng)
        @fastmath begin
            if Tr[1] > 10
                σ = min(max(σ*2.86f0*Accept[1]/Tr[1], 1/sqrt(step+1000)), 10)     # if the acceptance is higher (lower) than .35, sigma should increase (decreased) to explore more (less) the phase space to get closer to .35
            end
            if Tr[3] > 0
                a = 1/(Tr[3]+2)
                σϕ = min(max(σϕ*2.27f0*min(max(Accept[3]/Tr[3], a), 1-a), 1/sqrt(step + 1000)), 10)     # if the acceptance is higher (lower) than .44, sigma should increase (decreased) to explore more (less) the phase space to get closer to .44
            end         # I could have used clamp(a,b,c) instead of min(max(a,c),b) but it is less efficient (around 50% slower)
            if Tr[4] > 0
                a = 1/(Tr[4]+2)
                σθ = min(max(σθ*2.27f0*min(max(Accept[4]/Tr[4], a), 1-a), 1/sqrt(step+10000)), 2)     # if the acceptance is higher (lower) than .44, sigma should increase (decreased) to explore more (less) the phase space to get closer to .44
            end
            Acceptance += Accept
            Try += Tr
        end
    end
    return lat, σ, σϕ, σθ, Acceptance, Try, Cluster_update
end

function MH(N::Int64, T::Float32, L::Int, D::Float32, burn::Int, Skip)     # Sampler i.e. main function
    rng = Xoshiro(rand(1:10000))
    lat = Initial_lattice(L, rng)
    σ::Float32 = .25f0; σϕ::Float32 = .4f0; σθ::Float32 = .25f0
    Acceptance = zeros(Int64, 4); Try = zeros(Int64, 4)
    Cluster_update = 0
    Nmeasurement = div(N - burn, Skip)      # Number of measurement (as some sweeps are not measured to save some time)
    Energies = Array{Float32}(undef, Nmeasurement)
    Mag  = Array{Float32}(undef, Nmeasurement, 3)
    corr = zeros(Float32, length(RowCorr(lat, L)))
    
    for k=1:N
        lat, σ, σϕ, σθ, Acceptance, Try, Cluster_update = MHStep(k, lat, T, L, D, σ, σϕ, σθ, Acceptance, Try, Cluster_update, rng)
        if k > burn && mod(k, Skip) == 0            # to not measure every lattice sweeps
            m = div(k-burn, Skip)
            Energies[m] = Energy(lat, L, D)
            Mag[m, :] = mag(lat, L)
            corr += RowCorr(lat, L)
        end
    end
    Energies /= L^2
    corr /= Nmeasurement
    accept = Acceptance ./ Try
    @save "Data/[$L]_$N/$(L)_$(T)_$D.jld2" Energies Mx=Mag[:,1] My=Mag[:,2] Mz=Mag[:,3] corr accept
    println("$L \t $D \t $T \t $(round.(accept; digits=3)) \t and Magz2 : $(round(Mag_z2(getindex.(lat, 2), L); digits=3)) \t $(Cluster_update)")
    return Energies[end]
end

function RowCorr(Lattice::Array{SVector{2,Float32},2}, L::Int)   # Return the average row Correlation of the matrix (in a Correlation vector : first neigbour average correlation, second neigbour...)
    N=round(Int, L/2-1)
    corr = zeros(Float32, N)
    @inbounds @fastmath for i=1:N
        a=0f0
        for k=1:L
            k2 = (k + i <= L) ? k + i : k + i - L # like modulo, but more efficient (only if x < 2y for mod(x,y))
            for j=1:L
                a += Correlation(Lattice[j,k], Lattice[j,k2]) #cos(Lattice[1,j,k]-Lattice[1,j,mod(k+i-1,L)+1]) * cos(Lattice[2,j,k]-Lattice[2,j,mod(k+i-1,L)+1])
            end
        end
        corr[i] = a/L^2
    end
    return corr
end

function mag(Lattice::Array{SVector{2,Float32},2}, L::Int64)    # Compute the magnetization in x, y and z
    X = 0f0
    Y = 0f0
    Z = 0f0
    @inbounds @fastmath for k in 1:L, j in 1:L # faster as Julia stores arrays in column-major order
        lat_ij = Lattice[j,k]
        stheta, ctheta = sincos(lat_ij[2])
        sphi, cphi = sincos(lat_ij[1])
        x = cphi*stheta
        y = sphi*stheta
        z = ctheta
        X +=x
        Y +=y
        Z +=z
    end
    return [X, Y, Z]/L^2
end

function Mag_z2(ZLattice::Array{Float32,2}, L::Int64)   # Compute the sum(S_z^2), to check how how much the spin are align along the z-direction
    return sum(x -> cos(x)^2, ZLattice)/(L*L)
end

@inline @fastmath function IsingProba(Mz2::Float32)   # Probability of Ising-flip, just a choice, there is no perfect one
    return (1.9f0*Mz2 -.95f0)^2 *(sign(Mz2 - .5f0)+1)/2
end

@inline @fastmath function XYProba(Mz2::Float32)      # Probability of XY-flip (either phi or theta, independently): useful to converge faster for XY config as σϕ and σθ are tuned independently, to really have θ very close to π/2
    return (.9f0 - 9f0*Mz2)^2 *(sign(.1f0 -Mz2)+1)/2
end # just a choice, there is no perfect one

#    ----    get variables    ----    #

function parse_commandline()    # Parse command-line arguments
    s = ArgParseSettings()
    @add_arg_table! s begin
        # Required/Primary arguments
        "--n", "-n"
            help = "Number of sweeps"
            arg_type = Int64
            required = true

        "--l", "-l"
            help = "Lattice sizes (comma-separated, e.g., 8,12)"
            arg_type = String
            required = true
        
        "--t","-t"
            help = "Temperature (comma-separated, e.g., 0.63,0.67)"
            arg_type = String
            required = true
        
        "--d", "-d"
            help = "D values (comma-separated, e.g., -1,1)"
            arg_type = String
            required = true

        "--skip"
            help = "measure every \"skip\""
            arg_type = Int64
            default = 10
    end
    return parse_args(s)
end

function parse_csv(str::String, ::Type{T}) where T# Helper function to parse comma-separated values
    return [parse(T, x) for x in split(strip(str), ',')]
end


#   @inline     directly put the code of the function where the function is called, instead of really calling it    for single line function, no need "return" inside it
#   @fastmath   allows to reorder operation to save a bit of time (sometimes less precise because of floating point handling)
#   @inbounds   don't check if wanted component exist example     x=vector[n]     Doesn't check that the n-th component even exists, just trust the process