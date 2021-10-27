#
# in PKg mode, ]dev  Wick
using Revise

using Wick
using SymEngine
using JITFunc
# we add function to SSATape
ssatape=SSATape()
setInputArgs(ssatape,[:x,:y,:z])
r=evalTape(Basic("x+y"),ssatape)
evalTape(r,:r,ssatape)
#  we provide a shortcut to this function
r=ssatape(Basic("x+y"))
ssatape(r,:r)

# now, we can compile it as
func=SSACompiledFunc(ssatape,[:r])

x=rand(3);func(x...)[1]-(x[1]+x[2])

using JLD                       #  save the information

# function saveFunc(func::SSACompiledFunc,filename)
#     saveFunc(func.func,"$(filename)_func.bin")
#     save("$(filename)_info.jld","numberArray",func.numberArray,"computeArraySize",length(func.computeArray),"inputSize",func.inputSize,"outputSize",func.outputSize,"outputIndex",func.outputIndex)
# end

# load("$(filename)_info.jld")
filename="./ssacompiled_test"
saveSSAFunc(func,filename)
func=loadSSAFunc(filename)
func(rand(3)...)

#
r=ssatape(Basic(1.01))
convert(Basic,r,ssatape)

