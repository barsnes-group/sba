using Combinatorics
using LinearAlgebra
using Dates
include("preallocation.jl")
include("doptim.jl")

# Call this function to find an optimal allocation for a given set of
# sample sizes and batch sizes. Can take a long time if there are many
# samples, many subjects per sample, many batches, or any combination.
#
# samplesizes and batchsizes should be 1-dimensional column vectors
# 
# Returns the incidence matrix of an optimal allocation as a
# 2-dimensional Array of size (B, T) where B is the number of batches
# and T is the number of treatments, with the number of subjects of
# each treatment per batch.
function getOptimalAllocation(samplesizes, batchsizes, timing)
    topleft = makeTopLeft(samplesizes)
    bottomright = makeBottomRight(batchsizes)
    pre = preallocation!(samplesizes, batchsizes)
    optimalAllocation(zeros(Int, (0,length(samplesizes))), samplesizes, batchsizes, Combinatorics.combinations, topleft, bottomright, pre, timing)
end

# This is a helper function, do not call by hand.
function optimalAllocation(allocation, samplesizes, batchsizes, fun, topleft, bottomright, preallocation, timing)
    time0 = now()
    currentbatch = size(allocation, 1)+1
    allocated = dropdims(sum(allocation, dims=1), dims=1)
    subsleft = samplesizes .- allocated
    batchleft = length(batchsizes) - size(allocation, 1)
    if any(allocated .> samplesizes) || any(subsleft .> batchleft)
        return zeros(Int, (length(batchsizes), length(samplesizes)))
    elseif sum(batchsizes[currentbatch:length(batchsizes)]) == 0
        return vcat(allocation, zeros(Int, (batchleft, length(samplesizes))))+preallocation
    elseif size(allocation, 1) == length(batchsizes)
        return allocation+preallocation
    end
    deterministic = findall(subsleft .== batchleft)
    choiceset = filter(x -> x ∉ deterministic, findall(allocated .< samplesizes))
    choicen = batchsizes[currentbatch] - length(deterministic)
    batchiterator = fun(choiceset, choicen)
    iterlen = binomial(length(choiceset), choicen)
    if timing
        println(iterlen, " iterations to perform.")
    end
    bestdet = 0
    bestalloc = zeros(Int, (length(batchsizes), length(samplesizes)))
    for (idx, cmb) in enumerate(batchiterator)
        if timing
            if idx == 1
                println(now())
            else
                elapsed = now() - time0
                if (elapsed < Second(60))
                    println(round(idx-1 / iterlen, 1), "%\t", round(elapsed, Second(1)), " seconds.")
                else
                    println(round(idx-1 / iterlen, 1), "%\t", round(elapsed, Minute(1)), " minutes.")
                end
                estleft = elapsed/(idx-1) * (iterlen-idx+1)
                println("Estimated time left: ", Dates.format(estleft, "HH:MM"))
            end
        end
        newbatch = zeros(Int, (1, length(samplesizes)))
        newbatch[cmb] .= 1
        newbatch[deterministic] .= 1
        tempalloc = optimalAllocation(vcat(allocation, newbatch), samplesizes, batchsizes, fun, topleft, bottomright, preallocation, false)
        tempdet = doptim(tempalloc, topleft, bottomright)
        if tempdet > bestdet
            bestalloc = tempalloc
            bestdet = tempdet
        end
    end
    return bestalloc
end
