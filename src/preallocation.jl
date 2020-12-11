function preallocation(samplesizes::Array{<:Integer}, batchsizes::Array{<:Integer})
    return (preallocation!(copy(samplesizes),
                           copy(batchsizes)))
end

function preallocation!(samplesizes::Array{<:Integer}, batchsizes::Array{<:Integer})
    Nbatch = length(batchsizes)
    allocations = zeros(Int, (Nbatch,length(samplesizes)))
    for i in eachindex(samplesizes)
        nPerBatch = div(samplesizes[i], Nbatch)
        if nPerBatch > 0
            allocations[:, i] .+= nPerBatch
            samplesizes[i] = rem(samplesizes[i], Nbatch)
            batchsizes .-= nPerBatch
        end
    end
    allocations
end

