# implement the differentiation of the Mop to a Uop

# dev Wick

using Revise
using Wick
using SymEngine

op1=createOp("a_1",4)

op2=MultiOp(Basic("u"),UniOp("a_d_1",4))
op3=MultiOp(Basic("v"),UniOp("a_1",4))
op4=op2*op3

opDiff(op1,op1)

target=UniOp([1,3,4,5],4)
arg=UniOp([3,4],4)
target_=getOp(target)
arg_=getOp(arg)[1]
target_
arg_
opDiff(target_,4)

