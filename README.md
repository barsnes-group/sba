# Structured Stochastic Batch Allocation

Julia code for the allocation of a fixed cohort to batches.
Currently it is not a Module, but rather just a set of files.

First install [Julia](https://julialang.org) and download the repository.

## Optimal allocation
Run the following, either as a file or interactively.
The sample sizes and batch sizes are examples.
Take care that the sum of the batch sizes equals the sum of the sample sizes.

    include("optimalAllocation.jl")
    samplesizes = [5, 6, 7, 8, 9, 9]
    batchsizes = [8, 8, 7, 7, 7, 7]
    getOptimalAllocation(samplesizes, batchsizes)

This can take a very long time for larger cohorts.

## Quick random allocation
To run a quick random allocation of a larger cohort, when it takes too long to find the optimal allocation, run the following.
This runs both sba and the completely random binary allocation algorithms and picks the best out of a number of repeated runs.

    include("randomAllocation.jl")
    samplesizes = [5, 6, 7, 8, 9, 9]
    batchsizes = [8, 8, 7, 7, 7, 7]
    nruns = 1000
    getRandomAllocation(samplesizes, batchsizes, 1000)


## Quick allocation with a specific algorithm

### sba
To run a quick structured allocation of a larger cohort with sba, when it takes too long to find the optimal allocation, run the following.

    include("sba.jl")
    samplesizes = [5, 6, 7, 8, 9, 9]
    batchsizes = [8, 8, 7, 7, 7, 7]
    sba(samplesizes, batchsizes)

Given the speed and stochastic nature of the algorithm it is advised to run the algorithm multiple times and take the best one.

    include("doptim.jl")
    include("sba.jl")
    samplesizes = [5, 6, 7, 8, 9, 9]
    batchsizes = [8, 8, 7, 7, 7, 7]
    dopt = 0
    allocation = zeros(Int, (length(samplesizes), length(batchsizes)))
    topleft = makeTopLeft(samplesizes)
    bottomright = makeBottomRight(batchsizes)
    for i in 1:1000
        temp = hbmv(samplesizes, batchsizes)
        dtemp = doptim(temp, topleft, bottomright)
        if dopt < dtemp
            dopt = dtemp
            allocation = temp
        end
    end
    allocation

## Completely random (binary) allocation
To run a quick random allocation of a larger cohort, when it takes too long to find the optimal allocation, run the following.

    include("random.jl")
    samplesizes = [5, 6, 7, 8, 9, 9]
    batchsizes = [8, 8, 7, 7, 7, 7]
    randombinary(samplesizes, batchsizes)

Given the speed and stochastic nature of the algorithm it is advised to run the algorithm multiple times and take the best one.

    include("doptim.jl")
    include("random.jl")
    samplesizes = [5, 6, 7, 8, 9, 9]
    batchsizes = [8, 8, 7, 7, 7, 7]
    dopt = 0
    allocation = zeros(Int, (length(samplesizes), length(batchsizes)))
    topleft = makeTopLeft(samplesizes)
    bottomright = makeBottomRight(batchsizes)
    for i in 1:1000
        temp = randombinary(samplesizes, batchsizes)
        dtemp = doptim(temp, topleft, bottomright)
        if dopt < dtemp
            dopt = dtemp
            allocation = temp
        end
    end
    allocation
