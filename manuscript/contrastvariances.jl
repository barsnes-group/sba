using LinearAlgebra
using DelimitedFiles
using Random

include("../src/sba.jl")
using .SBA

include("random.jl")

# N is incidence matrix
function exampleincidence(optimal=true)
    if optimal
        [0 1 1 1 0; # batch 1
         1 0 0 1 1;
         0 1 1 0 1;
         1 1 1 0 0;
         1 0 0 1 1; # batch 5
         1 1 0 0 1;
         1 0 1 1 0;
         0 1 0 1 1;
         1 0 1 0 1;
         0 1 1 1 0; #batch 10
         0 0 1 1 1;
         0 1 0 1 1;
         0 0 1 0 1]
    else
        [1 1 1 0 0; # batch 1
         1 1 0 1 0;
         1 1 0 0 1;
         1 0 1 1 0;
         1 0 1 0 1; # batch 5
         1 0 0 1 1;
         0 1 1 1 0;
         0 1 1 0 1;
         0 1 0 1 1;
         0 0 1 1 1; #batch 10
         0 1 0 1 1;
         0 0 1 1 1;
         0 0 1 0 1]
    end
end

function exampleincidenceBailey()
    [1 1 1 1 0 0; 1 1 0 0 1 1; 0 0 1 1 1 1]
end

# this only works for a binary design!
function incidence2design(incidence)
    design = zeros(Int, (sum(incidence), size(incidence, 2)))
    i = 1
    for row=1:size(incidence, 1)
        for col=findall(incidence[row,:] .== 1)
            design[i,col] = 1
            i += 1
        end
    end
    design
end

function incidence2block(incidence)
    blocks = zeros(Int, (sum(incidence), sum(incidence)))
    counter = 1
    for bs=sum(incidence, dims=2)
        blocks[counter:(counter+bs-1),counter:(counter+bs-1)] .= 1
        counter += bs
    end
    blocks
end

function q(incidence)
    blocks = zeros(Float64, (sum(incidence), sum(incidence)))
    counter = 1
    for bs=sum(incidence, dims=2)
        blocks[counter:(counter+bs-1),counter:(counter+bs-1)] .= 1/bs
        counter += bs
    end
    LinearAlgebra.I(size(blocks,1)) - blocks
end


function l(incidence)
    incidence2design(incidence)' * q(incidence) * incidence2design(incidence)
end

function getContrastVariance(contrast, incidence)
    z = l(incidence) \ contrast
    z' * incidence2design(incidence)' * q(incidence) * incidence2design(incidence) * z
end

function allPairContrasts(incidence)
    Linv = pinv(l(incidence))
    rv = fill(NaN, (size(incidence,2), size(incidence,2)))
    for i=1:(size(incidence,2)-1)
        for j=(i+1):size(incidence,2)
            rv[i,j] = Linv[i,i] - Linv[i,j] - Linv[j,i] + Linv[j,j]
        end
    end
    rv
end

function writevartocsv(samplesizes, batchsizes, inner, outer, filename)
    towrite = fill("", (outer * 2*binomial(length(samplesizes), 2), 3))
    for fun in [rba, sba]
        funoffset = 0
        if fun==sba
            funoffset = binomial(length(samplesizes), 2) * outer
        end
        Random.seed!(1234)
        for run=1:outer
            pcs = allPairContrasts(fun(samplesizes, batchsizes, tracebreak=0, maxreps=inner))
            runoffset = (run-1) * binomial(length(samplesizes), 2)
            varoffset = 0
            for i=1:(length(samplesizes)-1)
                for j=(i+1):length(samplesizes)
                    varoffset += 1
                    towrite[runoffset + funoffset + varoffset, 1] = string(fun)
                    towrite[runoffset + funoffset + varoffset, 2] = string(i) * "-" * string(j)
                    towrite[runoffset + funoffset + varoffset, 3] = string(round(pcs[i,j], digits=5))
                end
            end
        end
    end
    open(filename, "w") do thefile
        write(thefile, "function,contrast,variance\n")
        writedlm(thefile, towrite, ",")
    end
    towrite 
end

function contrastvariances(inner=1000, outer=1000, postfix="contrastvar")
    writevartocsv([6,7,8,8,9], [3,3,3,3,3,3,3,3,3,3,3,3,2], inner, outer, "67889_3"*postfix*".csv") # C
end
