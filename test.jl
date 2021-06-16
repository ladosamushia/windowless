using FITSIO
using DataFrames
using QuadGK

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

"""
    Read eBOSS catalog from fits fits
"""
function read_fits_catalog(fitsfile)
    f = FITS(fitsfile)
    df = DataFrame(f[2])
    return df
end

function spherical_to_cartesian(ra, dec, z)
    ϕ = ra*π/180
    θ = dec*π/180 .+ π/2
    distance(z) = 2999/sqrt(0.3*(1+z)^3+0.7)
    r = zeros(length(z))
    for i in 1:length(r)
        r[i] = quadgk(distance, 0, z[i])[1]
    end
    x = r.*cos.(ϕ).*sin.(θ)
    y = r.*sin.(ϕ).*sin.(θ)
    z = r.*cos.(θ)
    return x, y, z
end

"""
    grid eBOSS data.
"""
function grid_survey(x, y, z, Ngrid)
    grid = zeros(Ngrid, Ngrid, Ngrid)
    xmin = minimum(x)
    ymin = minimum(y)
    zmin = minimum(z)
    Lx = maximum(x) - xmin
    Ly = maximum(y) - ymin
    Lz = maximum(z) - zmin
    L = maximum([Lx, Ly, Lz])
    for i in 1:length(x)
        ix = ceil(Int, (x[i] - xmin + 0.0001)/(L + 0.0001)*Ngrid)
        iy = ceil(Int, (y[i] - ymin + 0.0001)/(L + 0.0001)*Ngrid)
        iz = ceil(Int, (z[i] - zmin + 0.0001)/(L + 0.0001)*Ngrid)
        grid[ix, iy, iz] += 1
    end
    return grid, L/Ngrid
end