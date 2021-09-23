# in PKg mode, ]dev  Wick
using Revise

using Wick

# test DA and op
# L,N
da=DA(4,2)
typeof(op(da,1))
typeof(op(da,1.0))
# the order is l,t
op(da,"cr",[1,1])
op(da,"an",[2,1])
op1=op(da,"dm",[2,2,2,1])
typeof(op(op1,1.0))
typeof(op(op1,"1/2"))

# we now genreate the input for WIck thorem.
input=initInputMultiBandSpinSymmetric(da)
