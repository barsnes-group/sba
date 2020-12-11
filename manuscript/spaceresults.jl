using DelimitedFiles
using Random

include("../src/sba.jl")
using .SBA

# To get the results for timing and space
# optimalrun(6, "optimal4andon.csv")
# runspace(6, "mbs4andon.csv")


function optimalrun(maxsubs, filename, newseed=1234)
    Random.seed!(newseed)
    rv = zeros(Float64, (0, 5))
    for g1=3:maxsubs
        for g2=g1:maxsubs
            for g3=g2:maxsubs
                rv = vcat(rv, runbest(Int[g1, g2, g3]))
                for g4=g3:maxsubs
                    rv = vcat(rv, runbest(Int[g1, g2, g3, g4]))
                    for g5=g4:maxsubs
                        rv = vcat(rv, runbest(Int[g1, g2, g3, g4, g5]))
                        for g6=g5:maxsubs
                            rv = vcat(rv, runbest(Int[g1, g2, g3, g4, g5, g6]))
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
    for mbs=4:min(sum(samplesizes), 10)
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
                rv = vcat(rv, runbothspace(Int[g1, g2, g3]))
                for g4=g3:maxsubs
                    rv = vcat(rv, runbothspace(Int[g1, g2, g3, g4]))
                    for g5=g4:maxsubs
                        rv = vcat(rv, runbothspace(Int[g1, g2, g3, g4, g5]))
                        for g6=g5:maxsubs
                            rv = vcat(rv, runbothspace(Int[g1, g2, g3, g4, g5, g6]))
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
        naive = @timed SBA.getspace(samplesizes, batchsizes, true)
        space = @timed SBA.getspace(samplesizes, batchsizes, false)
        rv = vcat(rv, [string(samplesizes) mbs sum(samplesizes) naive[1] space[1] naive[2] space[2]])
    end
    rv
end



