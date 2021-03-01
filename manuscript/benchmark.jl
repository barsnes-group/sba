using DelimitedFiles
using Random

include("../src/sba.jl")
using .SBA

include("random.jl")

function getDeterminants(samplesizes, batchsizes, f, nRepeats)
    topleft = SBA.makeTopLeft(samplesizes)
    bottomright = SBA.makeBottomRight(batchsizes)
    determinants = fill(-1.0, nRepeats)
    for i in 1:nRepeats
        allocation = f(copy(samplesizes), copy(batchsizes))
        determinants[i] = SBA.dcrit(allocation, topleft, bottomright)
    end
    determinants
end

function writealltocsv(samplesizes, batchsizes, nRepeats, filename)
    towrite = zeros(Float64, (nRepeats, 0))
    for fun in [randombinary, SBA.singlerun]
        Random.seed!(1234)
        towrite = hcat(towrite, getDeterminants(samplesizes, batchsizes, fun, nRepeats))
    end
    open(filename, "w") do thefile
        write(thefile, "RBA,SBA\n")
        writedlm(thefile, towrite, ",")
    end
    towrite
end

function runall(inner=1000, outer=1; postfix="")
    writealltocsv(fill(6,5), fill(3,10), inner*outer, "5times6subsinbs3"*postfix*".csv") # A
    writealltocsv([5,6,7,8,9,9], [8,8,7,7,7,7], inner*outer, "blockprex"*postfix*".csv") # D
    writealltocsv([5,5,8,8,10,10,12,12,15,15], fill(5,20), inner*outer, "5.5.8.8.10.10.12.12.15.15bs5"*postfix*".csv") # E
    writealltocsv(fill(10,10), fill(5,20), inner*outer, "10times10subsinbs5"*postfix*".csv") # B
    writealltocsv([6,7,8,8,9], [3,3,3,3,3,3,3,3,3,3,3,3,2], inner*outer, "67889_3"*postfix*".csv") # C
end

function runmarathontime(inner=1000, outer=1000; postfix="")
    marathontime(fill(6,5), fill(3,10), inner, outer, "5times6subsinbs3"*postfix*".csv") # A
    marathontime([5,6,7,8,9,9], [8,8,7,7,7,7], inner, outer, "blockprex"*postfix*".csv") # D
    marathontime([5,5,8,8,10,10,12,12,15,15], fill(5,20), inner, outer, "5.5.8.8.10.10.12.12.15.15bs5"*postfix*".csv") # E
    marathontime(fill(10,10), fill(5,20), inner, outer, "10times10subsinbs5"*postfix*".csv") # B
    marathontime([6,7,8,8,9], [3,3,3,3,3,3,3,3,3,3,3,3,2], inner, outer, "67889_3"*postfix*".csv") # C
end

function marathontime(samplesizes, batchsizes, inner, outer, filename)
    topleft = SBA.makeTopLeft(samplesizes)
    bottomright = SBA.makeBottomRight(batchsizes)
    towrite = zeros(Float64, (outer, 0))
    for fun in [rba, sba]
        thisrundet = []
        thisruntime = []
        Random.seed!(1234)
        for i=1:outer
            allo = @timed fun(samplesizes, batchsizes, tracebreak=0, maxreps=inner)
            push!(thisrundet, SBA.dcrit(allo[1], topleft, bottomright))
            push!(thisruntime, allo[2])
        end
        towrite = hcat(towrite, thisrundet, thisruntime)
    end
    open(filename, "w") do thefile
        write(thefile, "RBA,RBAt,SBA,SBAt\n")
        writedlm(thefile, towrite, ",")
    end
    towrite    
end
