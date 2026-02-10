using BenchmarkTools


function test1(N)
    for i=1:N
        a = rand()
    end
    return a
end

function test2(N)
    a = rand(10^7)
    return a[end]
end

N = 10^7

@btime test1($N)
@btime test2($N)