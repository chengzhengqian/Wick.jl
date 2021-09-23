include("UniOp.jl")
include("MultiOp.jl")
export DA, op, initInputMultiBandSpinSymmetric
# we add some convenient interface to create operator
using Match

"""
a structure to encode the system information for a discrete action. L is the spatial size and N is the time size.
Notice we 
"""
struct DA
    L::Int
    N::Int
end

"""
show the L and N size
"""
function Base.show(io::IO,da::DA)
    print(io,"DA[L=$(da.L),N=$(da.N)]")
end

"""
a interface to create operator 
con acutally is not necessary
"""
function op(da::DA,str::String,pos::Array{Int}=Array{Int,1}(),c=Basic(1))
    @match str begin
        "cr" =>  cr(c,pos...,da)
        "an" =>  an(c,pos...,da)
        "dm" =>  dm(c,pos...,da)
    end    
end

function op(mop::MultiOp{T,F},c::Any) where {T,F}
    createMOp(c,mop)
end

function op(da::DA,c::F) where {T,F}
    cons(c,da)
end

"""
we use the time major scheme.
Previous, this function return both idx and size. We should break this into two functions
"""
function getIdx(l,t,da::DA)
    if(l<1 || l> da.L)
        error("l should be within 1 and $(da.L)")
    end
    if(t<1 || t> da.N)
        error("t should be within 1 and $(da.N)")
    end
    # s=da.L*da.N
    idx=(t-1)*da.L+l
    idx
end

"""
override the base function.
"""
function Base.size(da::DA)
    s=da.L*da.N
end



# creation operator
# c, coeffiicent, l, orbital index, t time index

function cons(c,da::DA)
    s=size(da)
    createMOp(c,createOp(Array{Int,1}(),s))
end

function cr(c,l,t,da::DA)
    idx=getIdx(l,t,da)
    s=size(da)
    createMOp(c,createOp([2*idx-1],s))
end

function an(c,l,t,da::DA)
    idx=getIdx(l,t,da)
    s=size(da)
    createMOp(c,createOp([2*idx],s))
end

# ad_1*a_2, the index for 1,2 are given by l1,t1,l2,t2
# notice the UniOp is timed ordered, so we need appropriate sign 
function dm(c,l1,t1,l2,t2,da::DA)
    idx1=getIdx(l1,t1,da)
    idx2=getIdx(l2,t2,da)
    s=size(da)
    if(idx1>idx2)
        c*=-1
    end    
    createMOp(c,createOp([2*idx1-1,2*idx2],s))
end

# now, we start to write the Wick theorem part.
# We directly convert a calculation to a more organized and general structure to record. (The general structure should be like llvm ir, which is independent of the given problem)

# the first step, we prepare the input for the Green function
"""
this is the default ordering the the spatial index for N_orbital 
1_up, 1_dn, ..      N_orbital_dn
1      ,   2    , .....,   N_spin_orbital
orbital_idx = 1, .., N_orbital
spin_idx = 1, 2 (up, dn)
"""
function defaultSpatialIndex(orbital_idx, spin_idx)
    (orbital_idx-1)*2+spin_idx
end

"""
we use a single symbol for each entry of the green function.
the green function is of format 
"""
function initInputMultiBandSpinSymmetric(da::DA;defaultSymbol="g0")
    input=Dict{Any,Basic}()
    if(da.L%2!=0)
        error("the spatial dimension of DA is $(L) which should be even number!")
    end
    n_orbital=trunc(Int,(da.L/2))
    s=size(da)
    for orb in 1:n_orbital
        for i in 1:da.N
            for j in 1:da.N
                i_up=(i-1)*da.L+defaultSpatialIndex(orb,1)
                i_dn=(i-1)*da.L+defaultSpatialIndex(orb,2)
                j_up=(j-1)*da.L+defaultSpatialIndex(orb,1)
                j_dn=(j-1)*da.L+defaultSpatialIndex(orb,2)
                sign=1
                if(j<i)
                    sign=-1
                end
                input[createOp([crIdx(i_up),anIdx(j_up)],s)]=Basic("$(defaultSymbol)_$(orb)_$(i)_$(j)")*sign
                input[createOp([crIdx(i_dn),anIdx(j_dn)],s)]=Basic("$(defaultSymbol)_$(orb)_$(i)_$(j)")*sign
            end
        end
    end
    input
end

