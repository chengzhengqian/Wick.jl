# UniOp, represent a single produce of creation operator, 
export createOp,UniOp,getOp
using Base
using BitIntegers
# a single operator in 2M bits
# 256 is sufficient for most cases
# thie type of T is match from 2^(2M)
# we now add more size, so we need to update SizeType
SizeType=UInt64

#  M is the size of spin orbitals, and 2*M is the number of total spin,
#  the operator presentes  a_d_1,a_1,a_d_2,a_2.., a_d_M,a_M
# the bit is form left to right

struct UniOp{T<:Unsigned}
    val:: T
    M::SizeType
end

# there is no coefficient inforamtion for UniOp
# find the  appropreate type for given s
# we add large size to support large system size
#  it seems that BitIntegers could allow arbitrary size, but for now, we just use the following version
function getTypeForM(s::Number)
    if(0<s<=4)
        UInt8
    elseif(4<s<=8)
        UInt16
    elseif(8<s<=16)
        UInt32
    elseif(16<s<=32)
        UInt64
    elseif(32<s<=64)
        UInt128
    elseif(64<s<=128)
        UInt256
    elseif(128<s<=256)
        UInt512
    elseif(256<s<=512)
        UInt1024
    else
        error("s should not be greater than 512")
    end
end

function createOp(val::Number,s::Number,shift::Number=0)
    s=convert(SizeType,s)
    val=convert(getTypeForM(s),val)
    UniOp(val<<shift,s)
end

function createOp(pos::Array{T,1},s::Number) where {T<:Number}
    s=convert(SizeType,s)
    val=convert(getTypeForM(s),0x0)
    tag=convert(getTypeForM(s),0x1)
    for i in pos
        val=val|(tag<<(i-1))
    end
    val=convert(getTypeForM(s),val)
    UniOp(val,s)
end

# get the operator in postion idx, 0 or 1

function getOp(k::T,M::SizeType,idx::Number) where T
    if(idx>2*M || idx<1)
        error("idx must not be larger than 2*M or smaller than 1")
    else
        tag=convert(getTypeForM(M),0x1)
        (k&(tag<<(idx-1)))>>(idx-1)
    end    
end

function getOp(op::UniOp{T},idx::Number) where T
    getOp(op.val,op.M,idx)
end

function getOp(uop::UniOp{T}) where T
    result=Array{Int,1}()
    for i in 1:2*uop.M
        if(getOp(uop,i)>0)
            append!(result,i)
        end
    end
    result
end

subscriptNumbers=["₀","₁","₂","₃","₄","₅","₆","₇","₈","₉"]
function convertNumToSubscript(num::Number)
    join([subscriptNumbers[parse(UInt8,c)+1]  for c in "$(num)"])
end

"""
map the index in space time to the index used in UniOp
for a the creation operator bit
"""
function crIdx(idx)
    2*idx-1
end

"""
map the index in space time to the index used in UniOp
for a the annihilation operator bit
"""
function anIdx(idx)
    2*idx
end



"""
we need to fix the print, see JuliaDiff.jl
"""
function Base.show(io::IO,x::UniOp{T}) where T
    if(x.val==0)
        print(io,"Î")
    else
        for i in 1:x.M
            if(getOp(x,crIdx(i))>0)
                print(io,"â†$(convertNumToSubscript(i))")
            end
            if(getOp(x,anIdx(i))>0)
                print(io,"â$(convertNumToSubscript(i))")
            end
        end
    end
end

# this is more convient way to construct the UniOp
function createOp(str::String,s::Number)
    tokens=split(str,"_")
    if(size(tokens)[1]==1)
        createOp(0x0,s)
    else
        type=0                  # 0 for creation, 1 for annilation
        if(size(tokens)[1]==3)
            pos=parse(SizeType,tokens[3])
            type=0
        else
            pos=parse(SizeType,tokens[2])
            type=1
        end
        if(type==0)
            idx=crIdx(pos)      # creation type
        else
            idx=anIdx(pos)
        end        
        if(idx>2*s)
            error("the $(str) has an invalid index for a system with size as $(s)")
        end        
        createOp(0x1,s,idx-1)
    end
end

# we define a simple interface
function UniOp(str::String,s::Number)
    createOp(str,s)
end

function UniOp(pos::Array{T,1},s::Number) where {T<:Number}
    createOp(pos,s)
end
