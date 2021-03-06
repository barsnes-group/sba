"""
    getlogspace(samplesizes, batchsizes)

Return an approximation of the number of possible allocations, log10.

# Examples
```jldoctest
julia> getlogspace(fill(6,5), fill(5,6))
0.0

julia> getlogspace(fill(6,5), fill(3,10))
9.0

julia> getlogspace(fill(10,10), fill(5,20))
45.62661027484934
```
"""
function getlogspace(samplesizes::Array{<:Integer}, batchsizes::Array{<:Integer})
    pa = preallocation!(copy(samplesizes), copy(batchsizes))
    nzg = sum(samplesizes .- sum(pa, dims=1)' .> 0) # nonzero groups
    if nzg == 0
        return 0.0
    end
    apb = sum(pa, dims=2)[1] # allocations per batch
    nzbs = batchsizes[findall(batchsizes .> apb)] .- apb
    pop!(nzbs) # remove last element, as that batch always has only one option

    space = 0
    for i=eachindex(nzbs)
        space += log(10, binomial(nzg, nzbs[i]))
    end
    space
end

function getspace(samplesizes::Array{<:Integer}, batchsizes::Array{<:Integer}, naive=true)
    getspace!(copy(samplesizes), copy(batchsizes), naive)
end

function getspace!(samplesizes::Array{<:Integer}, batchsizes::Array{<:Integer}, naive=true)
    preallocation!(samplesizes, batchsizes)
    if (sum(samplesizes) == 0)
        return 1
    else
        news = Integer[]
        newb = Integer[]
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

function naivespace(samplesizes::Array{<:Integer}, batchsizes::Array{<:Integer})
    if length(batchsizes) == 1
        return 1
    else
        return binomial(length(samplesizes), batchsizes[1]) *
            naivespace(samplesizes, batchsizes[2:end])
    end
end

function exhaustivespace(allocation::Array{<:Integer}, samplesizes::Array{<:Integer}, batchsizes::Array{<:Integer})
    exhaustivespace!(allocation, samplesizes, batchsizes)
end

function exhaustivespace!(allocation::Array{<:Integer}, samplesizes::Array{<:Integer}, batchsizes::Array{<:Integer})
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
