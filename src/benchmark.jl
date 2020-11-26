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
            tempdet = bestdet
            bestalloc = allocation
        end
    end
    bestalloc
end

function writealltocsv(samplesizes, batchsizes, nRepeats, filename)
    towrite = zeros(Float64, (nRepeats, 0))
    for fun in [randombinary, randomnonbinary, sba, sbatwo, sbathree]
        towrite = hcat(towrite, getDeterminants(samplesizes, batchsizes, fun, nRepeats))
    end
    open(filename, "w") do thefile
    write(thefile, "randBin,randNonBin,SBA,SBAtwo,SBAthree\n")
    writedlm(thefile, towrite, ",")
    end
    towrite
end

function runall()
    nruns=1000
    writealltocsv(fill(6,5), fill(3,10), nruns, "output/5times6subsinbs3.csv") # A
#    writealltocsv(fill(10,10), fill(5,20), nruns, "10times10subsinbs5.csv") # D
#    writealltocsv([4,4,4,3,3,3,2,2], [5,5,5,5,5], nruns, "binaryimbalance.csv") # B
    writealltocsv([5,6,7,8,9,9], [8,8,7,7,7,7], nruns, "output/blockprex.csv") # D
#    writealltocsv(fill(5,10), fill(5,10), nruns, "output/10times5subsinbs5.csv") # B
    writealltocsv(fill(10,10), fill(5,20), nruns, "output/10times10subsinbs5.csv") # B
    writealltocsv([6,7,8,8,9], [3,3,3,3,3,3,3,3,3,3,3,3,2], nruns, "output/67889_3.csv") # C
end

# number of groups
# 3 - 6
# number of subjects per group
# 3 - 10
# max batchsize
# 2 - 10

function betterrun(maxsubs, nruns, filename, newseed)
    Random.seed!(newseed)
    rv = zeros(Float64, (0, 4))
    for g1=3:maxsubs
        for g2=g1:maxsubs
            for g3=g2:maxsubs
                rv = vcat(rv, runthem([g1, g2, g3], nruns))
                for g4=g3:maxsubs
                    rv = vcat(rv, runthem([g1, g2, g3, g4], nruns))
                    for g5=g4:maxsubs
                        rv = vcat(rv, runthem([g1, g2, g3, g4, g5], nruns))
                        for g6=g5:maxsubs
                            rv = vcat(rv, runthem([g1, g2, g3, g4, g5, g6],nruns))
                            println([g1 g2 g3 g4 g5 g6])
                        end
                    end
                end
            end
        end
    end
    open(filename, "w") do thefile
        write(thefile, "samplesizes,batchsize,SBA,randBin\n")
        writedlm(thefile, rv, ",")
    end
    rv
end

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
    for mbs=3:min(sum(samplesizes), 10)
        # [samplesizes, batchsize, nsubs, optimalD, time]
        print("nsubs ", sum(samplesizes), " mbs ", mbs)
        batchsizes=maxbatchsizetobatchsizes(sum(samplesizes), mbs)
        optd = @timed getOptimalAllocation(copy(samplesizes), copy(batchsizes), timing=false, allocation=false)
        println(" time: ", optd[2])
        rv = vcat(rv, [string(samplesizes) mbs sum(samplesizes) optd[1] optd[2]])
    end
    rv
end

function runthem(samplesizes, nruns)
    rv = zeros(Float64, (0, 4))
    for mbs=2:min(sum(samplesizes), 10)
        # [samplesizes, batchsize, random, sba]
        batchsizes=maxbatchsizetobatchsizes(sum(samplesizes), mbs)
        rv = vcat(rv, [string(samplesizes) mbs maximum(getDeterminants(samplesizes, batchsizes, sbathree, nruns)) maximum(getDeterminants(samplesizes, batchsizes, randombinary, nruns))])
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
    for mbs=2:min(sum(samplesizes), 10)
        println(samplesizes, mbs)
        # [samplesizes, batchsize, nsubs, naive, space, ntime, stime]
        batchsizes=maxbatchsizetobatchsizes(sum(samplesizes), mbs)
        naive = @timed getspace(samplesizes, batchsizes, true)
        space = @timed getspace(samplesizes, batchsizes, false)
        rv = vcat(rv, [string(samplesizes) mbs sum(samplesizes) naive[1] space[1] naive[2] space[2]])
    end
    rv
end
