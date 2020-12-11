# Structured Stochastic Batch Allocation

Julia code for the allocation of a fixed cohort to batches.

First install [Julia](https://julialang.org) and download the repository.
Currently it is implemented as a module, but not added to the General Registry.

To run the algorithm

    include("../src/sba.jl")
    using .SBA

    samplesizes = [5, 6, 7, 8, 9, 9]
    batchsizes = [8, 8, 7, 7, 7, 7]
    sba(samplesizes, batchsizes)

If you only know the maximum batch size

    batchsizes = maxbatchsizetobatchsizes(samplesizes, 8)

In the (perhaps odd) situation where you only know the number of batches

    batchsizes = nrbatchtobatchsizes(samplesizes, 6)

To check how much computations it would (very) roughly take to exhaustively calculate the best allocation

    getlogspace(samplesizes, batchsizes)

Between 1e7 - 1e8 might take a couple of minutes, above 1e8 the computing time needed to exhaustively calculate the best allocation increases very rapidly.

If it seems feasible to exhaustively calculate the best allocation

    getOptimalAllocation(samplesizes, batchsizes)
