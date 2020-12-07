include("doptim.jl")
include("HBMV.jl")
include("random.jl")

function getRandomAllocation(samplesizes, batchsizes, nruns=1000)
    dopt = Inf
    allocation = zeros(Int, (length(samplesizes), length(batchsizes)))
    topleft = makeTopLeft(samplesizes)
    bottomright = makeBottomRight(batchsizes)
    for i in 1:1000
        temp = sba(samplesizes, batchsizes)
        dtemp = doptim(temp, topleft, bottomright)
        if dopt > dtemp
            dopt = dtemp
            allocation = temp
        end
        temp = randombinary(samplesizes, batchsizes)
        dtemp = doptim(temp, topleft, bottomright)
        if dopt > dtemp
            dopt = dtemp
            allocation = temp
        end
    end
    allocation
end

function onerun(samplesizes, batchsizes, nruns=1000)
    dopt = Inf
    allocation = zeros(Int, (length(samplesizes), length(batchsizes)))
    topleft = makeTopLeft(samplesizes)
    bottomright = makeBottomRight(batchsizes)
    for i in 1:1000
        temp = sba(samplesizes, batchsizes)
        dtemp = doptim(temp, topleft, bottomright)
        if dopt > dtemp
            dopt = dtemp
            allocation = temp
        end
        temp = randombinary(samplesizes, batchsizes)
        dtemp = doptim(temp, topleft, bottomright)
        if dopt > dtemp
            dopt = dtemp
            allocation = temp
        end
    end
    allocation
end
