using Combinatorics
#using LinearAlgebra
using Dates
include("preallocation.jl")

function getspace(samplesizes, batchsizes, naive=true)
    getspace!(copy(samplesizes), copy(batchsizes), naive)
end

function getspace!(samplesizes, batchsizes, naive=true)
    preallocation!(samplesizes, batchsizes)
    if (sum(samplesizes) == 0)
        return 1
    else
        news = []
        newb = []
        for s in samplesizes
            if s > 0
                push!(news, s)
            end
        end
        for b in batchsizes
            if b > 0
                push!(newb, b)
            end
        end
        if naive
            return naivespace(news, newb)
        else
            return exhaustivespace!(zeros(Int, (0, length(news))), news, newb)
        end
    end
end

function naivespace(samplesizes, batchsizes)
    if length(batchsizes) == 1
        return 1
    else
        return binomial(length(samplesizes), batchsizes[1]) *
            naivespace(samplesizes, batchsizes[2:end])
    end
end

function exhaustivespace!(allocation, samplesizes, batchsizes)
    time0 = now()
    currentbatch = size(allocation, 1)+1
    allocated = dropdims(sum(allocation, dims=1), dims=1)
    subsleft = samplesizes .- allocated
    batchleft = length(batchsizes) - size(allocation, 1)
    #println("batchsizes", batchsizes)
    #println("currentbatch", currentbatch)
    #println(batchleft)
    if any(allocated .> samplesizes) || any(subsleft .> batchleft)
        return 0
    elseif size(allocation, 1) == length(batchsizes)-1
        return 1
    end
    deterministic = findall(subsleft .== batchleft)
    choiceset = filter(x -> x ∉ deterministic, findall(allocated .< samplesizes))
    choicen = batchsizes[currentbatch] - length(deterministic)
    batchiterator = Combinatorics.combinations(choiceset, choicen)
    iterlen = binomial(length(choiceset), choicen)
    nallocs = 0
    for (idx, cmb) in enumerate(batchiterator)
        newbatch = zeros(Int, (1, length(samplesizes)))
        newbatch[cmb] .= 1
        newbatch[deterministic] .= 1
        nallocs += exhaustivespace!(vcat(allocation, newbatch), samplesizes, batchsizes)
    end
    nallocs
end
