function maxbatchsizetobatchsizes(N::Integer, mbs::Integer)::Array{Integer}
    if mbs < 1 || N < 1
        error("N and max batchsize should both be at least 1")
    end
    nbatch = Int(ceil((N/mbs)))
    rv = [div(N, nbatch) for i in 1:nbatch]
    restbatch = rem(N, nbatch)
    if restbatch > 0
        rv[1:restbatch] .+= 1
    end
    rv
end

function maxbatchsizetobatchsizes(N::Array{<:Integer}, mbs::Integer)::Array{Integer}
    maxbatchsizetobatchsizes(sum(N), mbs)
end

function nrbatchtobatchsizes(N::Integer, nbatch::Integer)::Array{Integer}
    if nbatch < 1 || N < 1
        error("N and max batchsize should both be at least 1")
    end
    rv = [div(N, nbatch) for i in 1:nbatch]
    restbatch = rem(N, nbatch)
    if restbatch > 0
        rv[1:restbatch] .+= 1
    end
    rv
end

function nrbatchtobatchsizes(N::Array{<:Integer}, nbatch)::Array{Integer}
    nrbatchtobatchsizes(sum(N), nbatch)
end
