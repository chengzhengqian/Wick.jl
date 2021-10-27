# in PKg mode, ]dev  Wick
using Revise

using Wick
using SymEngine
# test DA and op
ssatape=SSATape()
[addTape(i,ssatape) for i in [:x,:y,:z]]

r=evalTape(Basic("x+y+x*z"),ssatape)
evalTape(r,:r,ssatape)

input_size=length(ssatape.input)
args=[InputNode(i-1) for i in 1:input_size]
dict=initDefaultDiffDict(args,ssatape)
@time evalDiff(args,dict,ssatape)
dict
# there are redundant idensity (1.0) in expression. FIx in in convert .
ssatape
# evalDiff(Symbol("∂r∂x"),dict,ssatape)
addToOutput(Symbol("∂r∂x"),dict,ssatape)
addToOutput(Symbol("∂r∂y"),dict,ssatape)
addToOutput(Symbol("∂r∂z"),dict,ssatape)
addToOutput(Symbol("∂r1∂x"),dict,ssatape)
addToOutput(Symbol("∂r1∂y"),dict,ssatape)
addToOutput(Symbol("∂r1∂z"),dict,ssatape)
targets=[:r,:r1]
args=getInputArgs(ssatape)
# new API
addToOutput(targets,args,dict,ssatape)
r=evalTape(Basic("x*y*z"),ssatape)
evalTape(r,:r1,ssatape)

# we need to dictionary to hold the differential values
i1=InputNode(1)
c1=ComputeNode(1)

Derivative(InputNode(1),ComputeNode(1))==Derivative(InputNode(1),ComputeNode(1))
hash(Derivative(InputNode(1),ComputeNode(1)))
i=Derivative(InputNode(1),ComputeNode(1))

dict=Dict{Derivative,Node}()
for i in 1:length(ssatape.input)
    for j in 1:length(ssatape.input)
        diff=Derivative(InputNode(i-1),InputNode(j-1))
        value=1.0
        if(i!=j)
            value=0.0
        end        
        dict[diff]=evalTape(Basic(value),ssatape)
    end
end

current_compute_size=length(ssatape.compute)
input_size=length(ssatape.input)
ssaform=getfromA(ssatape.compute,ComputeNode(0))
arg=InputNode(0)
evalDiff(ssaform,InputNode(2),dict,ssatape)

node=ComputeNode(0)
args=[InputNode(i-1) for i in 1:input_size]
evalDiff(ComputeNode(0),args,dict,ssatape)
evalDiff(ComputeNode(1),args,dict,ssatape)

# chekc the diff using Zygote
using SymEngine

ssatape=SSATape()
args=[:x,:y,:z,:w]
args_nodes=setInputArgs(ssatape,args)
expr1="x+y*z-w*x"
expr2="-x+x*y*z-w*x"
expr3="-x+y*4.123-z*123"
r1=ssatape(Basic("$(expr1)"))
ssatape(r1,:r1)
r2=ssatape(Basic("$(expr2)"))
ssatape(r2,:r2)
r3=ssatape(Basic("$(expr3)"))
ssatape(r3,:r3)
ssatape

dict=initDefaultDiffDict(args_nodes,ssatape)
evalDiff(args_nodes,dict,ssatape)
values,syms=addToOutput([:r1,:r2,:r3],args,dict,ssatape)

outputs=[:r1,:r2,:r3,syms...]


func=SSACompiledFunc(ssatape,outputs)
test_input=rand(4)
result=func(test_input...)

args_str=join(args,",")
output_str=join([expr1,expr2,expr3],",")
func_julia=eval(Meta.parse("($(args_str))->[$(output_str)]"))
result_1=func_julia(test_input...)
result[1:3]-result_1
using Zygote
test_f=x->func_julia(x[1],x[2],x[3],x[4])
# notice the transpose, columb majore
result_2=reshape(transpose(Zygote.jacobian(test_f,test_input)[1]),:)
result_2_=result[4:end]
result_2-result_2_

saveSSAFunc(func,"./test_ssafunc")
func=loadSSAFunc("./test_ssafunc")
x,y,z,w=test_input
result_2[1]==1-w
result_2[2]==-1-w+y*z
result_2_[5]==-1-w+y*z

result_2_[2]
result_2_[3]
result_2[2]==z
result_2[3]==y
ssatape.output[:∂r1∂y]

func_interp=SSAInterpreter(ssatape,outputs)
result_3=func_interp(test_input...)
result_3[3:end]-result_2
# the interpreted version is correct, the problem is that output maybe a compute note or a array
