# Structured Stochastic Batch Allocation

Julia code for the allocation of a fixed cohort to batches.

First install [Julia](https://julialang.org) and download the repository.
Currently it is implemented as a module, but not added to the General Registry.

To run the algorithm

    include("src/sba.jl")
    using .SBA

    samplesizes = [5, 6, 7, 8, 9, 9]
    batchsizes = [8, 8, 7, 7, 7, 7]
    sba(samplesizes, batchsizes)

Without optional arguments, allocations are generated until there has been no improvement in the last 1000 tries.

To change the number of tries without improvement to 10 000:

    sba(samplesizes, batchsizes, tracebreak = 10 000)

The maximum number of tries can be set with:

    sba(samplesizes, batchsizes, maxreps = 10 000)

This will stop generating new allocations when there has been no improvement in the last 1000 tries, or when there have been 10 000 tries already, whichever comes first.

To run the algorithm until there has been no improvement in the last 1234 tries, while not allowing more than 2345 tries in total, and with seed 3456 for reproducibility:

    sba(samplesizes, batchsizes, tracebreak = 1234, maxreps = 2345, seed = 3456)

If you only know the maximum batch size

    batchsizes = maxbatchsizetobatchsizes(samplesizes, 8)

In the (perhaps odd) situation where you only know the number of batches

    batchsizes = nrbatchtobatchsizes(samplesizes, 6)

To check how much computations it would (very) roughly take to exhaustively calculate the best allocation

    getlogspace(samplesizes, batchsizes)

Between 1e7 - 1e8 might take a couple of minutes, above 1e8 the computing time needed to exhaustively calculate the best allocation increases very rapidly.

If it seems feasible to exhaustively calculate the best allocation

    getOptimalAllocation(samplesizes, batchsizes)
