using StatsBase

function rba(samplesizes::Array{<:Integer}, batchsizes::Array{<:Integer}; tracebreak=1000, maxreps=0, seed=nothing)
    tracebreak = tracebreak > 0 ? tracebreak : Base.Inf
    maxreps = maxreps > 0 ? maxreps : Base.Inf
    if seed != nothing
        Random.seed!(seed)
    end
    topleft = SBA.makeTopLeft(samplesizes)
    bottomright = SBA.makeBottomRight(batchsizes)
    bestdet = 0.0
    bestallo = zeros(Int, (length(batchsizes), length(samplesizes)))
    itrace = 0
    ireps = 0
    while ireps < maxreps && itrace < tracebreak
        ireps += 1
        itrace += 1
        newallocation = randombinary!(copy(samplesizes), copy(batchsizes))
        newdeterminant = SBA.dcrit(newallocation, topleft, bottomright)
        if newdeterminant > bestdet
            bestdet = newdeterminant
            bestallo = newallocation
            itrace = 0
        end
    end
    bestallo
end

function randombinary(samplesizes::Array{<:Integer}, batchsizes::Array{<:Integer})
    randombinary!(copy(samplesizes), copy(batchsizes))
end

function randombinary!(samplesizes::Array{<:Integer}, batchsizes::Array{<:Integer})
    prealloc = SBA.preallocation!(samplesizes, batchsizes)
    if any(samplesizes .>= length(batchsizes)) ||
        sum(samplesizes) != sum(batchsizes)
        return zeros(Int, (length(batchsizes), length(samplesizes)))
    end
    while true
        allocation = zeros(Int, (length(batchsizes), length(samplesizes)))
        allocated = zeros(Int, length(samplesizes))
        for batch in eachindex(batchsizes)
            groups = findall(allocated .< samplesizes)
            if length(groups) < batchsizes[batch]
                break
            end
            newsample = StatsBase.sample(groups, batchsizes[batch],
                                         replace=false, ordered=true)
            allocated[newsample] .+= 1
            allocation[batch, newsample] .= 1
            if batch == length(batchsizes)
                return prealloc+allocation
            end
        end
    end
    return prealloc
end

function randomnonbinary(samplesizes::Array{<:Integer}, batchsizes::Array{<:Integer})
    if sum(samplesizes) != sum(batchsizes)
        return zeros(Int, (length(batchsizes), length(samplesizes)))
    end
    while true
        allocation = zeros(Int, (length(batchsizes), length(samplesizes)))
        allocated = zeros(Int, length(samplesizes))
        left = samplesizes
        for batch in eachindex(batchsizes)
            groups = findall(left .> 0)
            chosen = Array{Int}(undef, 0)
            times = Array{Int}(undef, 0)
            while true
                newsample = StatsBase.sample(groups, batchsizes[batch],
                                             replace=true, ordered=true)
                newcount = Dict{Int, Int}()
                for i in newsample
                    newcount[i] = get(newcount, i, 0)+1
                end
                chosen = unique(keys(newcount))
                times = values(newcount)
                if !any(left[chosen] - collect(Iterators.rest(values(newcount))) .< 0)
                    break
                end
            end
            allocated[chosen] .+= times
            allocation[batch, chosen] .= times
            left[chosen] .-= times
            if batch == length(batchsizes)
                return allocation
            end
        end
    end
    return zeros(Int, (length(batchsizes), length(samplesizes)))
end
