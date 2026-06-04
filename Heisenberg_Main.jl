using Distributions, Base.Threads, Random, JLD2, StaticArrays, Base.Filesystem, ArgParse   # To have random distrib of spin, parallelize, interpolate, Random, save data in compact file, to use SVector which are optimized, to create folder
include("./Heisenberg_Main_functions.jl")

# Extract and convert arguments
parsed_args = parse_commandline()
const N = parsed_args["n"]
const L = parse_csv(parsed_args["l"], Int)
T_min, T_max, nT = Float32.(parse_csv(parsed_args["t"], Float32))
nT = Int64(nT)
const D = Float32.(parse_csv(parsed_args["d"], Float32))
const burn = Int(min(N/4,100000))     # Burning period
const swap = parsed_args["swap"]
const Skip = parsed_args["skip"]

if T_min == 0   # the if is just a way to cheat and use differently the parameters to have more T value around Tpic
        vect_T = T_max .+ [-.1, -.05, -.03, -.017, -.01, -.005, 0, .005, .01, .017, .03, .05, .1]
    if nT == 0 || nT == 1
        const T = round.(Float32.(vect_T);digits=3)
    else
        if swap == 80
            error("multiple swapping, but no first one was given")
        else
            const T = round.(Float32.(vect_T);digits=3)[aa:nT:end]
        end
    end
else
    const T = Temperatures(T_min, T_max, nT)
end

# Handle boolean flags (ArgParse returns Bool or nothing)
const PBC = get(parsed_args, "pbc", false) && !get(parsed_args, "no_pbc", false)
if !haskey(parsed_args, "pbc") && !haskey(parsed_args, "no_pbc")
    const PBC = true
end
const Save = get(parsed_args, "save", false) && !get(parsed_args, "no_save", false)
if !haskey(parsed_args, "save") && !haskey(parsed_args, "no_save")
    const Save = true
end
const SaveLattices = get(parsed_args, "save_lattices", false)


println("  --  $N sweeps  --  $T")
starttime = time()
if !isdir("Data/$(L)_$N")   # create the folder if it does not exist
    mkpath("Data/$(L)_$N")
end
for d in D
    for l in L
        E1_end = MH_parallel_tempering(N, T, L, l, d, burn, PBC, Save, SaveLattices, Skip, swap)
    end
end

t = round(Int, time()-starttime)    # record and save time
if Save == true
    open("Data/$(L)_$N/elapsed_time.txt", "a") do file
        write(file, "$t\n$D\n$T\n\n")
    end
end
println(t)