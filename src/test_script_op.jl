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
da=DA(4,2)
input,input_args=initInputMultiBandSpinSymmetric(da)
# we have  define the operator
op1=op(da,"dm",[1,1,1,1])*op(da,"dm",[1,2,1,2])

# we first initialize the ssatape
ssatape=SSATape()
setInputArgs(ssatape,input_args)
typeof(collect(keys(input))[1])
# we need to store the map of all uniop necessary to compute the given operator
# unlike the input, the symbolic expression correspoinds the  a node in ssatape or a numeric constance.
uniopMap=Dict{UniOp,Basic}()
for (k,v) in input
    # print("$(k)->$(v)\n")
    result=evalTape(v,ssatape)
    uniopMap[k]=convert(Basic,result,ssatape)
end


uop=UniOp((collect(keys(op1.val))[1]),op1.M)
evalWick(uop,uniopMap,ssatape)
evalWick(op1,uniopMap,ssatape)

# organize the code
# we now allow only evaluate the coefficient of op
da=DA(4,2)
input,input_args=initInputMultiBandSpinSymmetric(da)
u_args=[:u,:v]
ssatape=SSATape()
setInputArgs(ssatape,[input_args...,u_args...])
uniopMap=initUniOpMap(input,ssatape)
op1=Basic("u")*op(da,"dm",[1,1,1,1])*op(da,"dm",[1,2,1,2])
op2=Basic("v")*op(da,"dm",[2,1,2,1])*op(da,"dm",[2,2,2,2])
op3=Basic("u+v")*op(da,"dm",[3,1,3,1])*op(da,"dm",[3,2,3,2])
op4=simplify(op1*op2*op3)
mop=op4
evalWick(op4,uniopMap,ssatape)

r=evalWick(op2*op1*op3,uniopMap,ssatape)
mop=op2*op1*op3
val=collect(mop.val)[1][1]
uop=UniOp(val,mop.M)
getOp(uop)
mop=op1*op2+op1
mop_new=evalWick(mop,ssatape)
evalWick(mop_new,uniopMap,ssatape)
evalWick(mop,uniopMap,ssatape)

num_to_uop=Dict{Int,Set{UniOp}}()
for val in keys(mop.val)
    uop=UniOp(val,mop.M)
    num=length(getOp(uop))
    if(haskey(num_to_uop,num))
        push!(num_to_uop[num],uop)
    else
        num_to_uop[num]=Set{UniOp}([uop])
    end
end
uop.val

sort(collect(keys(num_to_uop)))
evalTape(r,:r1,ssatape)
r=evalWick(op2*op1,uniopMap,ssatape)
evalTape(r,:r2,ssatape)

outputSyms=[:r1,:r2]
func_interp=SSAInterpreter(ssatape,outputSyms)
input_test=rand(8)
@time func_interp(input_test...)

# we now test the compile module
asm=compile(ssatape)
using JITFunc
func_compiled=Func(asm, "test1",mode="nasm")
computeArray=initComputeArray(ssatape)
numberArray=initNumberArray(ssatape)
@time func_compiled(input_test,computeArray,numberArray)
getOutputArray(ssatape,outputSyms,computeArray)
# disassemble(func_compiled)
# to use the function, we need to generate
# test the new API
func=func_compiled=SSACompiledFunc(ssatape,outputSyms)
@time func_compiled(input_test...)


func_interp(input_test...)-getOutputArray(ssatape,outputSyms,computeArray)



# we now add differentiation function
#

#  the previous frequency changes the representation from Basic to Expr
# here, we try to only use Expr

mop1=MultiOp(1.0,UniOp([1,2],3))
mop2=MultiOp(2.0,UniOp([1,2],3))
typeof(mop1)
typeof(mop2)
mop1*mop2

# dump(:(1+2))
# a=Expr(:Number,1)

