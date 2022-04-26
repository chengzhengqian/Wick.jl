include("UniOp.jl")
include("MultiOp.jl")
export DA, op, initInputMultiBandSpinSymmetric,defaultSpatialIndex,initInput1dOneBandWithTranslationSymmetry, hopping, x_ops, initInput1dOneBandWithTranslationSymmetryAndZ, initInputGeneralN2, initInputGeneralN2WithZ
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

# we add some extra function to build operators, they are originally written for the two band model, and we move them here

"""
X operator for a given orbital and time step, i.e, the Hubbard oeprator for one band case
"""
function x_ops(da::DA,orb,t)
    idx_up=defaultSpatialIndex(orb,1)
    idx_dn=defaultSpatialIndex(orb,2)
    n_up=op(da,"dm",[idx_up,t,idx_up,t])
    n_dn=op(da,"dm",[idx_dn,t,idx_dn,t])
    p_up=1-n_up
    p_dn=1-n_dn
    Xe=p_up*p_dn
    Xup=n_up*p_dn
    Xdn=p_up*n_dn
    Xd=n_up*n_dn
    Xe,Xup,Xdn,Xd
end

"""
using MathExpr.jl as coefficient
"""
function x_ops(da::DA,orb,t,engine::ExprEngine)
    idx_up=defaultSpatialIndex(orb,1)
    idx_dn=defaultSpatialIndex(orb,2)
    n_up=op(da,"dm",[idx_up,t,idx_up,t],engine(1.0))
    n_dn=op(da,"dm",[idx_dn,t,idx_dn,t],engine(1.0))
    p_up=1-n_up
    p_dn=1-n_dn
    Xe=p_up*p_dn
    Xup=n_up*p_dn
    Xdn=p_up*n_dn
    Xd=n_up*n_dn
    Xe,Xup,Xdn,Xd
end

"""
a_i_spin1†a_j_spin2
1, spin up, 2 spin dn
hopping between two spin orbital, 
We update so we can take two different tiem step
we also add DA as parameters
"""
function hopping(da::DA,i,spin1,j,spin2,t1,t2)
    idx1=defaultSpatialIndex(i,spin1)
    idx2=defaultSpatialIndex(j,spin2)
    op(da,"dm",[idx1,t1,idx2,t2])
end

function hopping(da::DA,i,spin1,j,spin2,t1,t2,engine::ExprEngine)
    idx1=defaultSpatialIndex(i,spin1)
    idx2=defaultSpatialIndex(j,spin2)
    op(da,"dm",[idx1,t1,idx2,t2],engine(1.0))
end


"""
a_i_spin1†a_j_spin2
hopping between two spin orbital at a given tie step
"""
function hopping(da::DA,i,spin1,j,spin2,t)
    # idx1=defaultSpatialIndex(i,spin1)
    # idx2=defaultSpatialIndex(j,spin2)
    # op(da,"dm",[idx1,t,idx2,t])
    hopping(da,i,spin1,j,spin2,t,t)
end

function hopping(da::DA,i,spin1,j,spin2,t,engine::ExprEngine)
    # idx1=defaultSpatialIndex(i,spin1)
    # idx2=defaultSpatialIndex(j,spin2)
    # op(da,"dm",[idx1,t,idx2,t])
    hopping(da,i,spin1,j,spin2,t,t,engine)
end

include("./Input.jl")
