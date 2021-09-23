# in PKg mode, ]dev  Wick
using Revise

using Wick
a=BDict{Int,String}()

haskey(a.a_to_b,"12")
a=InputNode(1)
a=OutputNode(1)
a=ComputeNode(1)
c=BDict{Int,String}()
add(c,4=>"4")
getfromA(c,4)

# saform=SSAForm(ComputeNode(1),:+,[InputNode(1),InputNode(2)])
# ssaform1=SSAForm(ComputeNode(1),:-,[InputNode(1),InputNode(2)])
# ssaform2=SSAForm(ComputeNode(1),:-,[InputNode(1),InputNode(2)])

# ssaform1==ssaform2
# [InputNode(1)]==[InputNode(1)]
# (:-)==(:-)
# SSAForm(ComputeNode(1),:-,[InputNode(1),InputNode(2)])==SSAForm(ComputeNode(1),:-,[InputNode(1),InputNode(2)])

ssaform=SSAForm(:+,[InputNode(1),InputNode(2)])
SSAForm(:+,[InputNode(1),InputNode(2)])==SSAForm(:+,[InputNode(1),InputNode(2)])
a=Dict{SSAForm,Int}()
a[SSAForm(:+,[InputNode(1),InputNode(2)])]=2
a[SSAForm(:+,[InputNode(1),InputNode(2)])]=3
hash(SSAForm(:+,[InputNode(1),InputNode(2)]))
hash([InputNode(1),InputNode(3)])
hash([:-,InputNode(1),InputNode(3)])

# now, we write function to manipulate the ssatape
evalTape(:y,ssatape)
typeof(:(x+y))
dump(:(x+y*x))
expr=:(x+y*x)
expr.head
using Base.Meta
Meta.parse("x+y*x")
sym="I112"
ssaform=SSAForm(:+,[InputNode(1),InputNode(2)])
evalTape(ssaform,ssatape)
evalTape(:(x+y*x+x*x),ssatape)
evalTape(:I12,ssatape)
evalTape(:C12,ssatape)
evalTape(:O12,ssatape)
dump(:(1+2*x))
#
evalTape(1.2,ssatape)
ssatape=SSATape()
evalTape(:(x+y*z),ssatape)
evalTape(:(x+y*z+1.0*2.0),ssatape)
using SymEngine
expr=Basic("x+y+1.0*2.0")
evalTape(Meta.parse("$(expr)"),ssatape)
