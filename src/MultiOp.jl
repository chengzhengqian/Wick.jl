# represent a generic operator
#  stored as  f_i*O_i, where f is the coefficient and O_i is a UniOp
export MultiOp, simplify
using SymEngine                 # when F is Basic (i.e an symbolic expression)

struct MultiOp{T<:Unsigned,F}
    val:: Dict{T,F}
    M::SizeType
end


function Base.show(io::IO,x::MultiOp{T,F}) where {T,F}
    is_first=true
    if(length(x.val)==0)
        print(io,"0")
    else
        idx=0
        max_idx=10
        for (k,c) in x.val
            if(!is_first)
                print(io,"+")
            end
            if(idx>max_idx)
                print(io,"...")
                return
            end
            print(io,"(");print(io,c);print(io,")")
            print(io,UniOp(k,x.M))
            is_first=false
            idx+=1
        end
    end    
end

function createMOp(c::F,op::UniOp{T}) where {T,F}
    MultiOp(Dict(op.val=>c),op.M)
end


function createMOp(c::F,str::String,s::Number) where F
    op=createOp(str,s)
    createMOp(c,op)
end
# this creates a constance
function createMOp(c::F,s::Number)  where F
    op=createOp(0x0,s)
    createMOp(c,op)
end
# this creates a constance based on MultiOp
function createMOp(c::Any,mop::MultiOp{T,F})  where {T,F}
    op=createOp(0x0,mop.M)
    createMOp(convert(F,c),op)
end

# createMOp(1,p1)

function addMOpUOp!(mp::MultiOp{T,F},k::T,c::F) where {T,F}    
    if(haskey(mp.val,k))
        mp.val[k]+=c
        if(mp.val[k]==0)
            delete!(mp.val,k)
        end        
    else        
        mp.val[k]=c
    end
    mp
end

function addMOpUOp!(mp::MultiOp{T,F},sp::UniOp{T},c::F) where {T,F}
    addMOpUOp!(mp,sp.val,c)
end

function addMOpMOp(mp1::MultiOp{T,F},mp2::MultiOp{T,F}) where {T,F}
    result=MultiOp(Dict{T,F}(),mp1.M)
    for (k,c) in mp1.val
        addMOpUOp!(result,k,c)
    end    
    for (k,c) in mp2.val
        addMOpUOp!(result,k,c)
    end
    result
end

# this will add mp2 to mp1
function addMOpMOp!(mp1::MultiOp{T,F},mp2::MultiOp{T,F}) where {T,F}
    result=mp1
    for (k,c) in mp2.val
        addMOpUOp!(result,k,c)
    end
    result
end


# Now, we need to implement the multiplication
function getNum(sp::UniOp)
    sum([getOp(sp,i) for i in 1:2*sp.M])
end


#  the result is stored in MOp
"""
this compute the produce of two UniOp (k1,s) (k2,s), notice the non-trivial part is the result is not a UniOp, but generally a MultiOp. This follows from the ordering of how basis. i.e, a_1*a_d_1, and a_d_1*a_1 .
"""
function mulUOpUOp(k1::T,k2::T,s::Number,c::F) where {T,F}
    val=convert(T,0x0)
    s=convert(SizeType,s)
    result=MultiOp(Dict{T,F}(val=>c),s)
    size1=getNum(UniOp(k1,s))
    for i  in 1:s
        op1_cr=getOp(k1,s,crIdx(i))
        op1_an=getOp(k1,s,anIdx(i))
        op2_cr=getOp(k2,s,crIdx(i))
        op2_an=getOp(k2,s,anIdx(i))
        size1-=op1_cr+op1_an
        sign_=(-1)^((op2_cr+op2_an)*size1)
        ops=oneEntryMap[[op1_cr,op1_an,op2_cr,op2_an]]
        new_result=MultiOp(Dict{T,F}(),s)
        for op_ in ops
            addMOpMOp!(new_result, transformMOp(result,i,op_,sign_))
        end
        result=new_result
    end
    result
end


# i=1

# [op1_cr,op1_an,op2_cr,op2_an]
# [op_cr,op_an,sign]
# this should has at most two entries
oneEntryMap=Dict{Array{Int,1},Array{Array{Int,1},1}}(
    [0,0,0,0]=>[[0,0,1]],
    [0,0,0,1]=>[[0,1,1]],
    [0,0,1,0]=>[[1,0,1]],
    [0,1,0,0]=>[[0,1,1]],
    [1,0,0,0]=>[[1,0,1]],
    [0,0,1,1]=>[[1,1,1]],
    [1,1,0,0]=>[[1,1,1]],
    [1,0,1,0]=>[],
    [0,1,0,1]=>[],
    [1,0,0,1]=>[[1,1,1]],
    [0,1,1,0]=>[[1,1,-1],[0,0,1]],
    [1,1,1,0]=>[[1,0,1]],
    [1,1,0,1]=>[],
    [1,0,1,1]=>[],
    [0,1,1,1]=>[[0,1,1]],
    [1,1,1,1]=>[[1,1,1]],
)

# this is an auxillary function to implement the product
# n1,n2 cr an

# Be careful that the default data type of 0x is UInt8
"""
this set the bit at pos for all terms in mop, with a sign shit extraSign for all them
"""
function transformMOp(mop::MultiOp{T,F},pos::Number,data::Array{Int,1},extraSign=1) where {T,F}
    n1,n2,sign_=data
    idx_n1=pos*2-1
    idx_n2=pos*2
    result=MultiOp(Dict{T,F}(),mop.M)
    tag=convert(T,0x1)
    for (k,c) in mop.val
        if(n1==1)
            k=k|(tag<<(idx_n1-1))
        end
        if(n2==1)
            k=k|(tag<<(idx_n2-1))
        end
        c=c*sign_*extraSign
        addMOpUOp!(result,createOp(k,mop.M),c)
    end
    result
end

function mulMOpMOp(mp1::MultiOp{T,F},mp2::MultiOp{T,F}) where {T,F}
    if(mp1.M!=mp2.M)
        error("mp1 and mp2 must have the same size")
    else
        result=MultiOp(Dict{T,F}(),mp1.M)
        for (k1,c1) in mp1.val
            for (k2,c2) in mp2.val
                op=mulUOpUOp(k1,k2,mp1.M,c1*c2)
                addMOpMOp!(result,op)                
            end    
        end
        result
    end
end

#  this will update the given mp
# To make it consistent, we shoudl ensure n is also F type
function mulMOpNumber!(mp::MultiOp{T,F},n::F) where {T,F}
    for (k,c) in mp.val
        mp.val[k]=c*n
    end
    mp
end
# this will return a new mop
function mulMOpNumber(mp::MultiOp{T,F},n::F) where {T,F}
    result=MultiOp(Dict{T,F}(),mp.M)
    for (k,c) in mp.val
        result.val[k]=c*n
    end
    result
end

# we need to add the function to convert a number to MultiOp based on a givne mop, add to createMop

# add convient function to create MultiOp
function MultiOp(c::F,op::UniOp{T}) where {T,F}
    createMOp(c,op)
end

function MultiOp(c::F,str::String,s::Number) where F
    createMOp(c,str,s)
end
# this creates a constance
function MultiOp(c::F,s::Number)  where F
    createMOp(c,s)
end
# this creates a constance based on MultiOp
function MultiOp(c::Number,mop::MultiOp{T,F})  where {T,F}
    createMOp(c,mop)
end


# we now add the basic operator overload for MultiOp
# this is taken from Op.jl previous, Now we merge the relavent function to this file.
# one of the caveat in the previous code is that when we muliple a MultiOp{T,F} with a Number, we assume the number ::Number, but more accurately, we should assume it is F. i.e, we should deliberately enforce the code to use the same type in coefficient.
# using SymEngine
# typeof(Basic("1")*1)
# convert(Basic,1)
# convert(Basic,Basic("1"))
# maybe it is better to add convert in the code, we first change the Number to F and see

function Base.:*(mp1::MultiOp{T,F},mp2::MultiOp{T,F}) where {T,F}
    mulMOpMOp(mp1,mp2)
end

function Base.:*(mp1::MultiOp{T,F},c::F) where {T,F}
    mulMOpNumber(mp1,c)
end

function Base.:*(c::F,mp1::MultiOp{T,F}) where {T,F}
    mulMOpNumber(mp1,c)
end
# now, add the interface when the c is a generic Number, we first convert it to F, actually, it is Any, as string can also be converted into Basic
function Base.:*(mp1::MultiOp{T,F},c::Any) where {T,F}
    mp1*convert(F,c)
end

function Base.:*(c::Any,mp1::MultiOp{T,F}) where {T,F}
    mp1*convert(F,c)
end


function Base.:-(mp1::MultiOp{T,F}) where {T,F}
    mulMOpNumber(mp1,convert(F,-1))
end

function Base.:-(mp1::MultiOp{T,F},mp2::MultiOp{T,F}) where {T,F}
    addMOpMOp(mp1,mulMOpNumber(mp2,convert(F,-1)))
end

function Base.:-(mp::MultiOp{T,F},c::F) where {T,F}
    addMOpMOp(mp,createMOp(-c,mp))
end

function Base.:-(c::F,mp::MultiOp{T,F}) where {T,F}
    addMOpMOp(createMOp(c,mp),mulMOpNumber(mp,convert(F,-1)))
end

function Base.:-(mp::MultiOp{T,F},c::Any) where {T,F}
    mp-convert(F,c)
end

function Base.:-(c::Any,mp::MultiOp{T,F}) where {T,F}
    convert(F,c)-mp
end

function Base.:+(mp1::MultiOp{T,F},mp2::MultiOp{T,F}) where {T,F}
    addMOpMOp(mp1,mp2)
end

function Base.:+(mp::MultiOp{T,F},c::F) where {T,F}
    addMOpMOp(mp,createMOp(c,mp))
end

function Base.:+(c::F,mp::MultiOp{T,F}) where {T,F}
    addMOpMOp(mp,createMOp(c,mp))
end

function Base.:+(mp::MultiOp{T,F},c::Any) where {T,F}
    mp+convert(F,c)
end

function Base.:+(c::Any,mp::MultiOp{T,F}) where {T,F}
    mp+convert(F,c)
end

"""
by default, symbol engine does not expand the expression
this is only necessary when we the field is symbolic
"""
function simplify(mop::MultiOp{T,Basic}) where T
    for (k,c) in mop.val
        new_val=expand(c)
        if(new_val==0 || new_val==0.0)
            delete!(mop.val,k)
        else
            mop.val[k]=new_val
        end
    end
    mop
end

"""
simplify the term when using MathExpr
"""
function simplify(mop::MultiOp{T,SymExpr}) where T
    for (k,c) in mop.val
        new_val=MathExpr.simplify(c)
        if(is_zero(new_val))
            delete!(mop.val,k)
        else
            mop.val[k]=new_val
        end
    end
    mop
end

"""
When we use the opDiff, we require the argument is UniOp.
While the most code interface is for MultiOp. When MultiOp has only one entry, we could convert the MultiOp to UniOp, with the coefficient
"""
function UniOp(mop::MultiOp)
    s=mop.M
    if(length(mop.val)==1)
        k=collect(keys(mop.val))[1]
        UniOp(k,s),mop.val[k]
    else
        error("$(mop) should only have one entry")
    end    
end
