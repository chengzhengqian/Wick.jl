# to generate the input information for various models


"""
we use a single symbol for each entry of the green function.
the green function is of format.
Also, to facilitate the use of ssa tape, we return the order of the symbol
defaultSymbol="g0"
orb=1
n_orbital=2
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
    input, [Symbol("$(defaultSymbol)_$(orb)_$(i)_$(j)") for orb in 1:n_orbital for j in 1:da.N for i in 1:da.N]
end

"""
using our new library, MathExpr.jl
"""
function initInputMultiBandSpinSymmetric(da::DA,engine::ExprEngine;defaultSymbol="g0")
    input=Dict{Any,SymExpr}()
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
                input[createOp([crIdx(i_up),anIdx(j_up)],s)]=engine(Symbol("$(defaultSymbol)_$(orb)_$(i)_$(j)"))*sign
                input[createOp([crIdx(i_dn),anIdx(j_dn)],s)]=engine(Symbol("$(defaultSymbol)_$(orb)_$(i)_$(j)"))*sign
            end
        end
    end
    input, [Symbol("$(defaultSymbol)_$(orb)_$(i)_$(j)") for orb in 1:n_orbital for j in 1:da.N for i in 1:da.N]
end

"""
n_size=10
i=1
j=2
distanceFor1dChain(11,1,10)
"""
function distanceFor1dChain(i,j,n_size)
    min(abs(j-i)%n_size,(n_size-abs(j-i))%n_size)
end

"""
values=[ engine(:t11) engine(:t12) ; engine(:t21) engine(:t22)]
setInputForGivenTimeBlock(da,1,2,2,2,values)
values are the subblock for g for given spacial index described by orb1,spin1; orb2, spin2
"""
function setInputForGivenTimeBlock(da::DA,orb1,orb2,spin1,spin2,values)
    result=Dict{Any,Any}()
    s=size(da)    
    for t1 in 1:da.N
        for t2 in 1:da.N
            value_t1_t2=values[t1,t2]
            orb1_t1=(t1-1)*da.L+defaultSpatialIndex(orb1,spin1)
            orb2_t2=(t2-1)*da.L+defaultSpatialIndex(orb2,spin2)
            sign=1
            if(orb2_t2<orb1_t1)
                sign=-1
            end
            result[createOp([crIdx(orb1_t1),anIdx(orb2_t2)],s)]=value_t1_t2*sign
        end
    end
    result
end

"""
for N=2 only
input,input_args=initInput1dOneBandWithTranslationSymmetry(da,engine)
"""
function initInput1dOneBandWithTranslationSymmetry(da::DA,engine::ExprEngine;defaultSymbol="n0")
    input=Dict{Any,SymExpr}()
    if(da.N!=2)
        error("only works for N=2, but get $(da.N)")
    end
    n_orbital=trunc(Int,(da.L/2))
    for orb1 in 1:n_orbital
        for orb2 in 1:n_orbital
            dist=distanceFor1dChain(orb1,orb2,n_orbital)
            n_orb1_orb2_sym="$(defaultSymbol)_$(dist)"
            n_orb1_orb2=engine(Symbol(n_orb1_orb2_sym))
            if(dist==0)
                g_block=[ n_orb1_orb2 (1.0-n_orb1_orb2) ;
                          (-n_orb1_orb2) n_orb1_orb2]
            else
                g_block=[ n_orb1_orb2 (-n_orb1_orb2) ;
                          (-n_orb1_orb2) n_orb1_orb2]
            end
            merge!(input,setInputForGivenTimeBlock(da,orb1,orb2,1,1,g_block))
            merge!(input,setInputForGivenTimeBlock(da,orb1,orb2,2,2,g_block))
        end
    end
    max_dist=trunc(Int,floor(n_orbital/2))
    input, [Symbol("$(defaultSymbol)_$(dist)") for dist in 0:max_dist]
end

"""
this one assume we have a local self-energy that renormalized the density with z and main the density to n0_0, 
g_k={
    { nk,  z*(1-nk)},
    {-z*nk,z^2*(nk-n0)+n0}
}
So when transfrom the real space,
g_i_i={
    { n0, z*(1-n0)},
    {-z*n0, n0}
}
For off  diagonal part, for  i≠j
g_i_j={
    {nij,-z*nij},
    {-z*nij,z^2*nij}
}
"""
function initInput1dOneBandWithTranslationSymmetryAndZ(da::DA,engine::ExprEngine;defaultSymbol="n0",zSymbol="z")
    input=Dict{Any,SymExpr}()
    if(da.N!=2)
        error("only works for N=2, but get $(da.N)")
    end
    n_orbital=trunc(Int,(da.L/2))
    z=engine(Symbol(zSymbol))
    for orb1 in 1:n_orbital
        for orb2 in 1:n_orbital
            dist=distanceFor1dChain(orb1,orb2,n_orbital)
            n_orb1_orb2_sym="$(defaultSymbol)_$(dist)"
            n_orb1_orb2=engine(Symbol(n_orb1_orb2_sym))
            if(dist==0)
                # diagonal
                g_block=[ n_orb1_orb2  (z*(1.0-n_orb1_orb2)) ;
                          (-z*n_orb1_orb2) n_orb1_orb2]
            else
                # off diagonal
                g_block=[ n_orb1_orb2 (-z*n_orb1_orb2) ;
                          (-z*n_orb1_orb2) (z*z*n_orb1_orb2)]
            end
            merge!(input,setInputForGivenTimeBlock(da,orb1,orb2,1,1,g_block))
            merge!(input,setInputForGivenTimeBlock(da,orb1,orb2,2,2,g_block))
        end
    end
    max_dist=trunc(Int,floor(n_orbital/2))
    input,[ [Symbol("$(defaultSymbol)_$(dist)") for dist in 0:max_dist]...,Symbol(zSymbol)]
end

"""
to perform the general expansion, we should be able to specify the single density matrix arbitrarly.
Here, we does not impose any symmetry
for N=2 only
"""
function initInputGeneralN2(da::DA,engine::ExprEngine;defaultSymbol="n0")
    input=Dict{Any,SymExpr}()
    if(da.N!=2)
        error("only works for N=2, but get $(da.N)")
    end
    n_orbital=trunc(Int,(da.L/2))
    for orb1 in 1:n_orbital
        for orb2 in 1:n_orbital
            dist=distanceFor1dChain(orb1,orb2,n_orbital)
            n_orb1_orb2_sym="$(defaultSymbol)_$(orb1)_$(orb2)"
            n_orb1_orb2=engine(Symbol(n_orb1_orb2_sym))
            if(dist==0)
                g_block=[ n_orb1_orb2 (1.0-n_orb1_orb2) ;
                          (-n_orb1_orb2) n_orb1_orb2]
            else
                g_block=[ n_orb1_orb2 (-n_orb1_orb2) ;
                          (-n_orb1_orb2) n_orb1_orb2]
            end
            # spin up (1) and spin down (2)
            merge!(input,setInputForGivenTimeBlock(da,orb1,orb2,1,1,g_block))
            merge!(input,setInputForGivenTimeBlock(da,orb1,orb2,2,2,g_block))
        end
    end
    # max_dist=trunc(Int,floor(n_orbital/2))
    input, [Symbol("$(defaultSymbol)_$(orb1)_$(orb2)") for orb2 in 1:n_orbital for orb1 in 1:n_orbital] 
end


function initInputGeneralN2WithZ(da::DA,engine::ExprEngine;defaultSymbol="n0",zSymbol="z")
    input=Dict{Any,SymExpr}()
    if(da.N!=2)
        error("only works for N=2, but get $(da.N)")
    end
    n_orbital=trunc(Int,(da.L/2))
    z=engine(Symbol(zSymbol))
    for orb1 in 1:n_orbital
        for orb2 in 1:n_orbital
            dist=distanceFor1dChain(orb1,orb2,n_orbital)
            n_orb1_orb2_sym="$(defaultSymbol)_$(orb1)_$(orb2)"
            # n_orb1_orb2_sym="$(defaultSymbol)_$(dist)"
            n_orb1_orb2=engine(Symbol(n_orb1_orb2_sym))
            if(dist==0)
                # diagonal
                g_block=[ n_orb1_orb2  (z*(1.0-n_orb1_orb2)) ;
                          (-z*n_orb1_orb2) n_orb1_orb2]
            else
                # off diagonal
                g_block=[ n_orb1_orb2 (-z*n_orb1_orb2) ;
                          (-z*n_orb1_orb2) (z*z*n_orb1_orb2)]
            end
            merge!(input,setInputForGivenTimeBlock(da,orb1,orb2,1,1,g_block))
            merge!(input,setInputForGivenTimeBlock(da,orb1,orb2,2,2,g_block))
        end
    end
    # max_dist=trunc(Int,floor(n_orbital/2))
    input,[ [Symbol("$(defaultSymbol)_$(orb1)_$(orb2)") for orb2 in 1:n_orbital for orb1 in 1:n_orbital]..., Symbol(zSymbol) ]
    # input,[ [Symbol("$(defaultSymbol)_$(dist)") for dist in 0:max_dist]...,Symbol(zSymbol)]
end

# reshape([(i,j) for j in 1:2 for i in 1:2],2,2)
