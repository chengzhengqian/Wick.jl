# in PKg mode, ]dev  Wick
using Revise

using Wick
using SymEngine

mop1=MultiOp(1,UniOp("a_d_1",4))
mop2=MultiOp(1,"a_d_1",4)
mop3=MultiOp(2,4)
mop4=MultiOp(4,mop3)


op1=MultiOp(Basic("u"),UniOp("a_d_1",4))
op2=MultiOp(Basic("v"),UniOp("a_1",4))
op3=op1*op2
op4=op2*op1
op3+op4
# we add overload for either F or Number, (now, we just use Any)
# convert(Basic,"1/2")
op5=op4*"1/2"                   #  now this code also works
op5-op3
op5-1.0
-1.0-op5
op5+1.0
1.0+op5
op5*1.0
1.0*op5-1*op5
simplify(op1+op2+op3)

