"""
    Create index of nonzero cells.
"""
function inon0(gridr)
    Ngrid = size(gridr)[1]
    # number of nonzero gridcells
    N = sum(gridr.!=0)
    # keep x, y, z indexes of nonzero gridcells
    index = zeros(Int, 3, N)
    counter = 1
    for ix in 1:Ngrid, iy in 1:Ngrid, iz in 1:Ngrid
        if gridr[ix, iy, iz] != 0
            index[:, counter] = [ix iy iz]
            counter += 1
        end
    end
    return index
end

"""
    distance between two sets of indeces.
"""
function distance(ix1, iy1, iz1, ix2, iy2, iz2)
    ceil(Int, sqrt((ix1 - ix2)^2 + (iy1 - iy2)^2 + (iz1 - iz2)^2))
end

"""
    loop over pairs of gridcells
"""
function pair_loop(index, pairs)
    Ngrid = size(index)[2]
    for i in 1:Ngrid 
        println(i)
        for j in i:Ngrid
            id = distance(index[1,i], index[2,i], index[3,i], index[1,j], index[2,j], index[3,j])
            if id == 0
                id += 1
            end
            pairs[id] += 1
        end
    end
    return pairs
end