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
opDiff(target_,3)
opDiff(target_,[1,3])
opDiff(target_,[3,1,2])

target=MultiOp(Basic("u*v"),UniOp([1,3,4],4))+MultiOp(Basic("u"),UniOp([1,3,4,5],4))
arg1=UniOp([1],4)
arg2=UniOp([4],4)
simplify(opDiff(target,arg))
opDiff(target,arg1,arg1)
opDiff(target,arg2,arg1,arg1)

da=DA(4,2)
mop=op(da,"cr",[1,2])
UniOp(mop)

marg=MultiOp(-1,arg)
mop=simplify(opDiff(target,4*marg))*4

typeof(mop)
mop.val[0x00]==0.0
