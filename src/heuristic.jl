function getStartPair(lambda::Array{<:Integer}, groups::Array{<:Integer})
    minimum=lambda[groups[1],groups[2]]
    minarr=[groups[1] groups[2]]
    for i=1:(length(groups)-1), j=(i+1):length(groups)
        if lambda[groups[i],groups[j]] == minimum
            minarr = vcat(minarr, [groups[i] groups[j]])
        elseif lambda[groups[i],groups[j]] < minimum
            minimum=lambda[groups[i],groups[j]]
            minarr=[groups[i] groups[j]]
        end
    end
    minarr[sample(1:size(minarr, 1), 1)[1], :]
end

function getStartSingle(lambda::Array{<:Integer}, groups::Array{<:Integer})
    groups[sample(findall(lambda[diagind(lambda)[groups]] .== minimum(lambda[diagind(lambda)[groups]])), 1)[1]]
end

"""
    sba(samplesizes, batchsizes; tracebreak=1000, maxreps=0, seed=nothing)

Return a batch allocation for the given sample and batch sizes.  To
keep trying allocations until no improvement is found for x tries, set
`tracebreak`.  To limit the total number of tries, set `maxreps`.
Setting either to anything less than 1 will set the given value to
infinity.

See also: maxbatchsizetobatchsizes, nrbatchtobatchsizes

# Examples
```jldoctest
julia> sba([6,6,6,6,6], [5,5,5,5,5,5])
6×5 Array{Int64,2}:
 1  1  1  1  1
 1  1  1  1  1
 1  1  1  1  1
 1  1  1  1  1
 1  1  1  1  1
 1  1  1  1  1

julia> sba(fill(6,5), fill(3,10), seed=123)
10×5 Array{Int64,2}:
 1  0  0  1  1
 1  1  1  0  0
 0  1  1  1  0
 0  0  1  1  1
 1  1  0  0  1
 0  1  1  0  1
 1  0  1  1  0
 0  1  0  1  1
 1  0  1  0  1
 1  1  0  1  0

julia> sba(fill(6,5), fill(3,10), tracebreak=10, seed=123)
10×5 Array{Int64,2}:
 1  1  0  0  1
 1  0  1  1  0
 0  1  1  0  1
 1  0  0  1  1
 0  1  1  1  0
 1  0  1  0  1
 1  1  0  1  0
 0  0  1  1  1
 0  1  0  1  1
 1  1  1  0  0

julia> sba(fill(6,6), fill(4,9), tracebreak=0, maxreps=20, seed=123)
9×6 Array{Int64,2}:
 1  0  1  1  0  1
 1  1  1  0  1  0
 0  1  0  1  1  1
 1  0  0  1  1  1
 0  1  1  1  0  1
 1  1  1  1  0  0
 0  1  1  0  1  1
 1  1  0  1  1  0
 1  0  1  0  1  1
```
"""
function sba(samplesizes::Array{<:Integer}, batchsizes::Array{<:Integer}; tracebreak=1000, maxreps=0, seed=nothing)
    if tracebreak <= 0 && maxreps <= 0
        throw(DomainError((tracebreak, maxreps), "either tracebreak or maxreps must be strictly positive"))
    end
    tracebreak = tracebreak > 0 ? tracebreak : Base.Inf
    maxreps = maxreps > 0 ? maxreps : Base.Inf
    if seed != nothing
        Random.seed!(seed)
    end
    topleft = makeTopLeft(samplesizes)
    bottomright = makeBottomRight(batchsizes)
    bestdet = 0.0
    bestallo = zeros(Int, (length(batchsizes), length(samplesizes)))
    itrace = 0
    ireps = 0
    while ireps < maxreps && itrace < tracebreak
        ireps += 1
        itrace += 1
        # println(ireps, "\t", itrace)
        newallocation = singlerun!(copy(samplesizes), copy(batchsizes))
        newdeterminant = dcrit(newallocation, topleft, bottomright)
        if newdeterminant > bestdet
            bestdet = newdeterminant
            bestallo = newallocation
            itrace = 0
        end
    end
    bestallo
end

function singlerun(samplesizes::Array{<:Integer}, batchsizes::Array{<:Integer})
    singlerun!(copy(samplesizes), copy(batchsizes))
end

function singlerun!(samplesizes::Array{<:Integer}, batchsizes::Array{<:Integer})
    # Initial allocation
    preallo = preallocation!(samplesizes, batchsizes)
    allocations = zeros(Int, (length(batchsizes), length(samplesizes)))
    lambda = zeros(Int, (length(samplesizes), length(samplesizes)))
    # Rest allocation
    for bsi in eachindex(batchsizes)
        groups = findall(samplesizes .> 0)
        # deterministic allocations
        chosen = findall(samplesizes .== length(batchsizes) - (bsi-1))
        neverChosen = findall(dropdims(sum(lambda, dims=1), dims=1) .== 0)
        neverChosen = filter(x -> x ∈ groups, neverChosen)
        neverChosen = filter(x -> x ∉ chosen, neverChosen)
        if length(neverChosen) > 0
            if length(neverChosen) > batchsizes[bsi] - length(chosen)
                chosen = vcat(chosen, StatsBase.sample(
                    neverChosen, batchsizes[bsi] - length(chosen),
                    replace=false, ordered=true))
            else
                chosen = vcat(chosen, neverChosen)
                while length(chosen) < batchsizes[bsi]
                    colsums = dropdims(sum(view(lambda, chosen,
                                                filter(x -> x ∉ chosen, groups)), dims=1), dims=1)
                    nextchoice = sample(findall(colsums .== minimum(colsums)), 1)[1]
                    chosen = vcat(chosen, filter(x -> x ∉ chosen, groups)[nextchoice])
                end
            end
        elseif length(chosen) < batchsizes[bsi] # if length(neverChosen)
            if length(chosen) == 0 # pick first pair of treatments
                if batchsizes[bsi] == 1
                    chosen = vcat(chosen, getStartSingle(lambda, groups))
                else
                    chosen = getStartPair(lambda, groups)
                end
            end
            while length(chosen) < batchsizes[bsi]
                colsums = dropdims(sum(view(lambda, chosen, filter(x -> x ∉ chosen, groups)), dims=1), dims=1)
                nextchoice = sample(findall(colsums .== minimum(colsums)), 1)[1]
                chosen = vcat(chosen, filter(x -> x ∉ chosen, groups)[nextchoice])
            end
        end # if length(neverChosen) else
        allocations[bsi, chosen] .+= 1
        for rowi=chosen, coli=chosen
            lambda[rowi,coli] += 1
        end
        samplesizes[chosen] .-= 1
    end # for bsi in eachindex(batchsizes)
    return(allocations+preallo)
end

