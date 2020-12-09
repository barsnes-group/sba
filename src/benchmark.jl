using DelimitedFiles
using Random
include("doptim.jl")
include("batches.jl")
include("random.jl")
include("sba.jl")
include("searchspace.jl")

function getDeterminants(samplesizes, batchsizes, f, nRepeats)
    topleft = makeTopLeft(samplesizes)
    bottomright = makeBottomRight(batchsizes)
    determinants = fill(-1.0, nRepeats)
    for i in 1:nRepeats
        allocation = f(copy(samplesizes), copy(batchsizes))
        determinants[i] = doptim(allocation, topleft, bottomright)
    end
    determinants
end

function getBest(samplesizes, batchsizes, f, nRepeats)
    topleft = makeTopLeft(samplesizes)
    bottomright = makeBottomRight(batchsizes)
    bestdet = 0.0
    bestalloc = zeros(Int, (length(batchsizes), length(samplesizes)))
    for i in 1:nRepeats
        allocation = f(copy(samplesizes), copy(batchsizes))
        tempdet = doptim(allocation, topleft, bottomright)
        if tempdet > bestdet
            bestdet = tempdet
            bestalloc = allocation
        end
    end
    bestalloc
end

function writealltocsv(samplesizes, batchsizes, nRepeats, filename)
    towrite = zeros(Float64, (nRepeats, 0))
    for fun in [randombinary, sba]
        Random.seed!(1234)
        towrite = hcat(towrite, getDeterminants(samplesizes, batchsizes, fun, nRepeats))
    end
    open(filename, "w") do thefile
        write(thefile, "randBin,SBA\n")
        writedlm(thefile, towrite, ",")
    end
    towrite
end

function writemarathon(samplesizes, batchsizes, outer, inner, filename)
    funcs = [randombinary, sba]
    towrite = zeros(Float64, (inner*outer, length(funcs)))
    for (idx, fun) in enumerate(funcs)
        Random.seed!(1234)
        towrite[:, idx] = getDeterminants(samplesizes, batchsizes, fun, outer*inner)
    end
    open(filename, "w") do thefile
        write(thefile, "randBin,SBAthree\n")
        writedlm(thefile, towrite, ",")
    end
    towrite
end

function runall(inner=1000, outer=1, postfix="")
    writealltocsv(fill(6,5), fill(3,10), inner*outer, "5times6subsinbs3"*postfix*".csv") # A
#    writealltocsv([4,4,4,3,3,3,2,2], [5,5,5,5,5], nruns, "binaryimbalance.csv") # B
    writealltocsv([5,6,7,8,9,9], [8,8,7,7,7,7], inner*outer, "blockprex"*postfix*".csv") # D
    writealltocsv([5,5,8,8,10,10,12,12,15,15], fill(5,20), inner*outer, "5.5.8.8.10.10.12.12.15.15bs5"*postfix*".csv") # E
    writealltocsv(fill(10,10), fill(5,20), inner*outer, "10times10subsinbs5"*postfix*".csv") # B
    writealltocsv([6,7,8,8,9], [3,3,3,3,3,3,3,3,3,3,3,3,2], inner*outer, "67889_3"*postfix*".csv") # C
end

function runmarathon(inner=1000, outer=1, postfix="")
    marathons(fill(6,5), fill(3,10), inner, outer, "5times6subsinbs3"*postfix*".csv") # A
#    marathons([4,4,4,3,3,3,2,2], [5,5,5,5,5], nruns, "binaryimbalance.csv") # B
    marathons([5,6,7,8,9,9], [8,8,7,7,7,7], inner, outer, "blockprex"*postfix*".csv") # D
    marathons([5,5,8,8,10,10,12,12,15,15], fill(5,20), inner*outer, "5.5.8.8.10.10.12.12.15.15bs5"*postfix*".csv") # E
    marathons(fill(10,10), fill(5,20), inner, outer, "10times10subsinbs5"*postfix*".csv") # B
    marathons([6,7,8,8,9], [3,3,3,3,3,3,3,3,3,3,3,3,2], inner, outer, "67889_3"*postfix*".csv") # C
end

function runmarathontime(inner=1000, outer=1, postfix="")
    marathontime(fill(6,5), fill(3,10), inner, outer, "5times6subsinbs3"*postfix*".csv") # A
#    marathons([4,4,4,3,3,3,2,2], [5,5,5,5,5], nruns, "binaryimbalance.csv") # B
    marathontime([5,6,7,8,9,9], [8,8,7,7,7,7], inner, outer, "blockprex"*postfix*".csv") # D
    marathontime([5,5,8,8,10,10,12,12,15,15], fill(5,20), inner*outer, "5.5.8.8.10.10.12.12.15.15bs5"*postfix*".csv") # E
    marathontime(fill(10,10), fill(5,20), inner, outer, "10times10subsinbs5"*postfix*".csv") # B
    marathontime([6,7,8,8,9], [3,3,3,3,3,3,3,3,3,3,3,3,2], inner, outer, "67889_3"*postfix*".csv") # C
end

function marathontime(samplesizes)
    towrite = zeros(Float64, (outer, 0))
    for fun in [randombinary, sba]
        thisrundet = []
        thisruntime = []
        Random.seed!(1234)
        for i=1:outer
            allo = @timed getAllocation(samplesizes, batchsizes, inner, fun, false)
            push!(thisrundet, allo[1])
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

function marathons(samplesizes, batchsizes, inner, outer, filename)
    towrite = zeros(Float64, (outer, 0))
    for fun in [randombinary, sba]
        thisrun = []
        Random.seed!(1234)
        for i=1:outer
            push!(thisrun, getAllocation(samplesizes, batchsizes, inner, fun, false))
        end
        towrite = hcat(towrite, thisrun)
    end
    open(filename, "w") do thefile
        write(thefile, "randBin,SBA\n")
        writedlm(thefile, towrite, ",")
    end
    towrite    
end

function getAllocation(samplesizes, batchsizes, nreps=1000, fun=sba, allocation=true, seed=nothing)
    if seed != nothing
        Random.seed!(seed)
    end
    topleft = makeTopLeft(samplesizes)
    bottomright = makeBottomRight(batchsizes)
    bestdet = 0.0
    bestallo = zeros(Int, (length(batchsizes), length(samplesizes)))
    for i in 1:nreps
        newallocation = fun(copy(samplesizes), copy(batchsizes))
        newdeterminant = doptim(newallocation, topleft, bottomright)
        if newdeterminant > bestdet
            bestdet = newdeterminant
            bestallo = newallocation
        end
    end
    if allocation
        return bestallo
    else 
        return bestdet
    end
end

# number of groups
# 3 - 6
# number of subjects per group
# 3 - 10
# max batchsize
# 2 - 10

function optimalrun(maxsubs, filename, newseed)
    Random.seed!(newseed)
    rv = zeros(Float64, (0, 5))
    for g1=3:maxsubs
        for g2=g1:maxsubs
            for g3=g2:maxsubs
                rv = vcat(rv, runbest([g1, g2, g3]))
                for g4=g3:maxsubs
                    rv = vcat(rv, runbest([g1, g2, g3, g4]))
                    for g5=g4:maxsubs
                        rv = vcat(rv, runbest([g1, g2, g3, g4, g5]))
                        for g6=g5:maxsubs
                            rv = vcat(rv, runbest([g1, g2, g3, g4, g5, g6]))
                            println([g1 g2 g3 g4 g5 g6])
                        end
                    end
                end
            end
        end
    end
    open(filename, "w") do thefile
        write(thefile, "samplesizes,batchsize,nsubs,optimal,time\n")
        writedlm(thefile, rv, ",")
    end
    rv
end

function runbest(samplesizes)
    rv = zeros(Float64, (0, 5))
    for mbs=min(sum(samplesizes), 10)
        # [samplesizes, batchsize, nsubs, optimalD, time]
        print("nsubs ", sum(samplesizes), " mbs ", mbs)
        batchsizes=maxbatchsizetobatchsizes(sum(samplesizes), mbs)
        optd = @timed getOptimalAllocation(copy(samplesizes), copy(batchsizes), timing=false, allocation=false)
        println(" time: ", optd[2])
        rv = vcat(rv, [string(samplesizes) mbs sum(samplesizes) optd[1] optd[2]])
    end
    rv
end



function runspace(maxsubs, filename)
    rv = zeros(Float64, (0, 7))
    for g1=3:maxsubs
        for g2=g1:maxsubs
            for g3=g2:maxsubs
                rv = vcat(rv, runbothspace([g1, g2, g3]))
                for g4=g3:maxsubs
                    rv = vcat(rv, runbothspace([g1, g2, g3, g4]))
                    for g5=g4:maxsubs
                        rv = vcat(rv, runbothspace([g1, g2, g3, g4, g5]))
                        for g6=g5:maxsubs
                            rv = vcat(rv, runbothspace([g1, g2, g3, g4, g5, g6]))
                            println([g1 g2 g3 g4 g5 g6])
                        end
                    end
                end
            end
        end
    end
    open(filename, "w") do thefile
        write(thefile, "samplesizes,batchsize,nsubs,naive,space,ntime,stime\n")
        writedlm(thefile, rv, ",")
    end
    rv
end

function runbothspace(samplesizes)
    rv = zeros(Float64, (0, 7))
    for mbs=4:min(sum(samplesizes), 10)
        println(samplesizes, mbs)
        # [samplesizes, batchsize, nsubs, naive, space, ntime, stime]
        batchsizes=maxbatchsizetobatchsizes(sum(samplesizes), mbs)
        naive = @timed getspace(samplesizes, batchsizes, true)
        space = @timed getspace(samplesizes, batchsizes, false)
        rv = vcat(rv, [string(samplesizes) mbs sum(samplesizes) naive[1] space[1] naive[2] space[2]])
    end
    rv
end

function allnaive(maxsubs, filename)
    rv = zeros(Float64, (0, 6))
    for g1=3:maxsubs
        for g2=g1:maxsubs
            for g3=g2:maxsubs
                rv = vcat(rv, runnaivespace([g1, g2, g3]))
                for g4=g3:maxsubs
                    rv = vcat(rv, runnaivespace([g1, g2, g3, g4]))
                    for g5=g4:maxsubs
                        rv = vcat(rv, runnaivespace([g1, g2, g3, g4, g5]))
                        for g6=g5:maxsubs
                            rv = vcat(rv, runnaivespace([g1, g2, g3, g4, g5, g6]))
                            println([g1 g2 g3 g4 g5 g6])
                        end
                    end
                end
            end
        end
    end
    open(filename, "w") do thefile
        write(thefile, "samplesizes,batchsize,nsubs,ngroups,nbatch,naive\n")
        writedlm(thefile, rv, ",")
    end
    rv
end

function runnaivespace(samplesizes)
    rv = zeros(Float64, (0, 6))
    for mbs=2:min(sum(samplesizes), 10)
        # [samplesizes, batchsize, nsubs, ngroups, nbatch, naive]
        batchsizes=maxbatchsizetobatchsizes(sum(samplesizes), mbs)
        rv = vcat(rv, [string(samplesizes) mbs sum(samplesizes) length(samplesizes) length(batchsizes) getlogspace(samplesizes, batchsizes, 10)])
    end
    rv
end








##############################

function marathonvars(inner=1000, outer=1000, postfix="contrastvar")
    writevartocsv(fill(6,5), fill(3,10), inner, outer, "5times6subsinbs3"*postfix*".csv") # A
    #    writealltocsv([4,4,4,3,3,3,2,2], [5,5,5,5,5], nruns, "binaryimbalance.csv") # B old
    # writealltocsv([5,6,7,8,9,9], [8,8,7,7,7,7], inner*outer, "output/blockprex"*postfix*".csv") # D
    # writealltocsv(fill(10,10), fill(5,20), inner*outer, "output/10times10subsinbs5"*postfix*".csv") # B
    writevartocsv([6,7,8,8,9], [3,3,3,3,3,3,3,3,3,3,3,3,2], inner, outer, "67889_3"*postfix*".csv") # C
end

function writevartocsv(samplesizes, batchsizes, inner, outer, filename)
    towrite = fill("", (outer * 2*binomial(length(samplesizes), 2), 3))
    for fun in [randombinary, sba]
        funoffset = 0
        if fun==sba
            funoffset = binomial(length(samplesizes), 2) * outer
        end
        Random.seed!(1234)
        for run=1:outer
            pcs = allPairContrasts(getAllocation(samplesizes, batchsizes, inner, fun, true))
            runoffset = (run-1) * binomial(length(samplesizes), 2)
            varoffset = 0
            for i=1:(length(samplesizes)-1)
                for j=(i+1):length(samplesizes)
                    varoffset += 1
                    towrite[runoffset + funoffset + varoffset, 1] = string(fun)
                    towrite[runoffset + funoffset + varoffset, 2] = string(i) * "-" * string(j)
                    towrite[runoffset + funoffset + varoffset, 3] = string(round(pcs[i,j], digits=5))
                end
            end
        end
    end
    open(filename, "w") do thefile
        write(thefile, "function,contrast,variance\n")
        writedlm(thefile, towrite, ",")
    end
    towrite 
end
