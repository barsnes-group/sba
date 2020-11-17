function preallocation(samplesizes,
                       batchsizes)
    return (preallocation!(copy(samplesizes),
                           copy(batchsizes)))
end

function preallocation!(samplesizes, batchsizes)
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
