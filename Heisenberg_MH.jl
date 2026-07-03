using Distributions, Base.Threads, Random, JLD2, StaticArrays, Base.Filesystem, ArgParse   # To have random distrib of spin, parallelize, interpolate, Random, save data in compact file, to use SVector which are optimized, to create folder
include("./Heisenberg_MH_functions.jl")

# Extract and convert arguments
parsed_args = parse_commandline()
const N = parsed_args["n"]
const L = parse_csv(parsed_args["l"], Int)
const T = Float32.(parse_csv(parsed_args["t"], Float32))
const D = Float32.(parse_csv(parsed_args["d"], Float32))
const burn = Int(min(N/4,100000))     # Burning period
const Skip = parsed_args["skip"]


println("  --  $N sweeps  --  $T")
starttime = time()
if !isdir("Data/$(L)_$N")   # create the folder if it does not exist
    mkpath("Data/$(L)_$N")
end
for d in D
    for l in L
        for t in T
            E1_end = MH(N, t, l, d, burn, Skip)
        end
    end
end

t = round(Int, time()-starttime)    # record and save time
open("Data/$(L)_$N/elapsed_time.txt", "a") do file
    write(file, "$t\n$D\n$T\n\n")
end
println(t)