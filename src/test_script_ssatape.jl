# in PKg mode, ]dev  Wick
using Revise

using Wick
using SymEngine
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

# treat power and minux
dump(:(-x+1))
dump(:(x-1))
dump(:(x^2))
2 isa Number
2==2.0
evalTape(:(x^2),ssatape)
evalTape(:(x^3),ssatape)
r=evalTape(:(-x+1),ssatape)
evalTape(r,:r1,ssatape)
r=evalTape(Basic("0"),ssatape)
s=convert(Basic,r,ssatape)
typeof(s)
evalTape(s,ssatape)

# we update the outpu form, now test it.
using SymEngine
ssatape=SSATape()
r=evalTape(Basic("x+1+y*x+z*y"),ssatape)
# add result to output. (now, we don't specify the output order)
evalTape(r,:r1,ssatape)
outputArray=[:r1]
inputArray=rand(3)
interpreter=SSAInterpreter(ssatape,outputArray)

interpret(InputNode(0),interpreter,inputArray)
interpret(ComputeNode(0),interpreter,inputArray)
interpret(NumberNode(0),interpreter,inputArray)
interpreter(0.0,0.0,0.0)


# now, we write the intepreter for the SSATape
# we test a more complex example
expr1="x+y+z^2-w+1*(2-x+y*w)"
expr2="x*4-y+z^2-w+1*(2-x+y*w)+(-w)*x"
ssatape=SSATape()
r=evalTape(Basic(expr1),ssatape)
# add result to output. (now, we don't specify the output order)
evalTape(r,:r1,ssatape)
r=evalTape(Basic(expr2),ssatape)
# add result to output. (now, we don't specify the output order)
evalTape(r,:r2,ssatape)
outputArray=[:r1,:r2]
interpreter=SSAInterpreter(ssatape,outputArray)
inputArgs=getInputArgs(interpreter)
test_input=rand(length(inputArgs))
compiled_func=eval(Meta.parse("($(join(inputArgs,",")))->[$(expr1),$(expr2)]"))
interpreter(test_input...)-compiled_func(test_input...)
@time interpreter(test_input...)
@time compiled_func(test_input...)
# we can write the machine code generation after we finish the Wick theorm.
r=evalTape(Basic("0"),ssatape)
convert(Basic,r,ssatape)
