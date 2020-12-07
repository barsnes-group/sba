using LinearAlgebra

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

