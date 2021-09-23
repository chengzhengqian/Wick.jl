# in PKg mode, ]dev  Wick
using Revise

using Wick

# for UniOp.jl
# for annilation operator
op1=createOp("a_1",5)
# for creation operator (d for dagger, anything thing seems OK from the code)
op1=createOp("a_d_1",5)
op1=createOp("I",5)
# we add new interface

op1=UniOp("a_1",2)
op2=UniOp("a_2",2)
op3=UniOp("a_d_2",2)
op4=UniOp("a_3",2)               #  should thow a error
op5=UniOp("a_d_3",2)               #  should thow a error

# using Pkg
# Pkg.add("BitIntegers")
# using BitIntegers
# UInt512(0x1)
# UInt1024(0x1)<<1000
# c=UInt128(0x1)
# supertypes(UInt128)
op=UniOp("a_d_4",5)
op=UniOp("a_d_430",500)

UniOp([1,2,3],5)
UniOp([2,3],5)
