function maxbatchsizetobatchsizes(N, mbs)
    if mbs < 1 || N < 1
        error("N and max batchsize should both be at least 1")
    end
    nbatch = Int(ceil((N/mbs)))
    rv = [div(N, nbatch) for i in 1:nbatch]
    restbatch = rem(N, nbatch)
    if restbatch > 0
        rv[1:restbatch] .+= 1
    end
    return(rv)
end

function nrbatchtobatchsizes(N, nbatch)
    if nbatch < 1 || N < 1
        error("N and max batchsize should both be at least 1")
    end
    rv = [div(N, nbatch) for i in 1:nbatch]
    restbatch = rem(N, nbatch)
    if restbatch > 0
        rv[1:restbatch] .+= 1
    end
    return(rv)
end