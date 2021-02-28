module SBA

using StatsBase
using Random
using LinearAlgebra
using Dates
using Combinatorics

include("preallocation.jl")
include("searchspace.jl")
include("optimalAllocation.jl")
include("dcrit.jl")
include("batches.jl")
include("heuristic.jl")

export sba, getlogspace, getOptimalAllocation, maxbatchsizetobatchsizes, nrbatchtobatchsizes

end

