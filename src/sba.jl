using StatsBase
using LinearAlgebra
include("preallocation.jl")

function getStartPair(lambda, groups)
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

function getStartSingle(lambda, groups)
    groups[sample(findall(lambda[diagind(lambda)[groups]] .== minimum(lambda[diagind(lambda)[groups]])), 1)[1]]
end

function sba(samplesizes, batchsizes)
    sba(copy(samplesizes), copy(batchsizes))
end

function sba!(samplesizes, batchsizes)
    # Initial allocation
    allocations = preallocation!(samplesizes, batchsizes)
    lambda = fixLambda(allocations)
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
    return(allocations)
end
