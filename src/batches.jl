"""
    maxbatchsizetobatchsizes(N, mbs)

Return a one-dimensional Array with batch sizes from the maximum
allowed batch size and either a one-dimensional Array with sample
sizes or the total sample size. Batches are made as equal to each
other in size as possible.

See also: nrbatchtobatchsizes

# Examples
```jldoctest
julia> maxbatchsizetobatchsizes(25, 4)
7-element Array{Integer,1}:
 4
 4
 4
 4
 3
 3
 3

julia> maxbatchsizetobatchsizes([2,4,6,12], 4)
6-element Array{Integer,1}:
 4
 4
 4
 4
 4
 4
```
"""
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

"""
    nrbatchtobatchsizes(N, nbatch)

Return a one-dimensional Array with batch sizes from the requested
number of batches and either a one-dimensional Array with sample sizes
or the total sample size. Batches are made as equal to each other in
size as possible.

See also: maxbatchsizetobatchsizes

# Examples
```jldoctest
julia> nrbatchtobatchsizes(25, 4)
4-element Array{Integer,1}:
 7
 6
 6
 6

julia> nrbatchtobatchsizes([2,4,6,12], 4)
4-element Array{Integer,1}:
 6
 6
 6
 6
```
"""
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
