# a simplified version of SSA, non condition right now

export InputNode, OutputNode, ComputeNode, BDict, getfromA, getfromB,hasA, hasB, add,SSAForm, SSATape, evalTape

abstract type Node end

struct InputNode <: Node
    idx::Int
end

const ReservedInputNodeSymbol="I"
"""
input node
"""
function Base.show(io::IO,input::InputNode)
    print(io,"$(ReservedInputNodeSymbol)$(input.idx)")
end

struct OutputNode <: Node
    idx::Int
end

const ReservedOutputNodeSymbol="O"
"""
output node
"""
function Base.show(io::IO,output::OutputNode)
    print(io,"$(ReservedOutputNodeSymbol)$(output.idx)")
end

struct ComputeNode <: Node
    idx::Int
end

const ReservedComputeNodeSymbol="C"
"""
output node
"""
function Base.show(io::IO,compute::ComputeNode)
    print(io,"$(ReservedComputeNodeSymbol)$(compute.idx)")
end


struct NumberNode <: Node
    idx::Int
end

const ReservedNumberNodeSymbol="N"
"""
output node
"""
function Base.show(io::IO,number::NumberNode)
    print(io,"$(ReservedNumberNodeSymbol)$(number.idx)")
end

struct BDict{A,B}
    a_to_b::Dict{A,B}
    b_to_a::Dict{B,A}
end

function BDict{A,B}() where {A,B}
    BDict{A,B}(Dict{A,B}(),Dict{B,A}())
end

function Base.length(bdict::BDict{A,B}) where {A,B}
    length(bdict.a_to_b)
end


function Base.show(io::IO,bdict::BDict{A,B}) where {A,B}
    print(io,"$(A)<=>$(B):$(bdict.a_to_b)")
end


function getfromA(bdict::BDict{A,B},a::A) where {A,B}
    bdict.a_to_b[a]
end

function getfromB(bdict::BDict{A,B},b::B) where {A,B}
    bdict.b_to_a[b]
end

function hasA(bdict::BDict{A,B},a::A) where {A,B}
    haskey(bdict.a_to_b,a)
end

function hasB(bdict::BDict{A,B},b::B) where {A,B}
    haskey(bdict.b_to_a,b)
end

"""
the BDict should strictly be an one-to-one map
"""
function add(bdict::BDict{A,B},m::Pair{A,B}) where {A,B}
    a,b=m
    if(hasA(bdict,a) || hasB(bdict,b))
        error("the BDict should be an one-to-one map ")
    end
    bdict.a_to_b[a]=b
    bdict.b_to_a[b]=a
end

"""
just the operation plus the arguments
"""
struct SSAForm
    op::Symbol
    args::Array{Node,1}
end

function Base.show(io::IO,ssaform::SSAForm)
    print(io,"($(join(ssaform.args,ssaform.op)))")
end

function Base.:(==)(ssa1::SSAForm,ssa2::SSAForm)
    (ssa1.op==ssa2.op)&&(ssa1.args==ssa2.args)
end

function Base.isequal(ssa1::SSAForm,ssa2::SSAForm)
    (ssa1.op==ssa2.op)&&(ssa1.args==ssa2.args)
end

"""
this is necessary to identify the same ssa in the Dict.
"""
function Base.hash(ssa::SSAForm)
    hash([ssa.op,ssa.args...])
end


"""
We use BDict to store the information
"""
struct SSATape    
    input_name::BDict{InputNode,Symbol} #  track the order of the input args
    compute::BDict{ComputeNode,SSAForm}
    output::Dict{OutputNode,Node} # each output node strictly correpodes a compute node, but a compute node may correspond zero or multiple output node (or maybe just input node, though that does not make too much sense here)
    output_name::BDict{OutputNode,Symbol} #  track the order of the output args
    number::BDict{NumberNode,Float64}     #  we assume all number are Float64
end

function SSATape()
    SSATape(BDict{InputNode,Symbol}(),BDict{ComputeNode,SSAForm}(),
            Dict{OutputNode,ComputeNode}(),BDict{OutputNode,Symbol}(),
            BDict{NumberNode,Float64}())
end

"""
return the InputNode for the symbol.
If the symbol does not exists, create a new one
We follow the convection that Node idx start from 0;
This means that when we generate the assembly code, we can directly use this index. (though this is not consistent with the julia's default convention)
As we going to use SymEngine to perform some simplification. So if an symbol start with I, OR C, indicates that it is an input nodes or compute nodes.
We assume that if it start with I,C,O, the index makes sense, and we does not check right now.
"""
function evalTape(sym::Symbol,tape::SSATape)
    sym_str="$(sym)"
    if(startswith(sym_str,ReservedInputNodeSymbol))
        return InputNode(parse(Int,sym_str[2:end]))
    elseif(startswith(sym_str,ReservedOutputNodeSymbol))
        return OutputNode(parse(Int,sym_str[2:end]))
    elseif(startswith(sym_str,ReservedComputeNodeSymbol))
        return ComputeNode(parse(Int,sym_str[2:end]))
    elseif(startswith(sym_str,ReservedNumberNodeSymbol))
        return NumberNode(parse(Int,sym_str[2:end]))
    else
        input=tape.input_name
        if(!hasB(input,sym))
            newNode=InputNode(length(input))
            add(input,newNode=>sym)
        end    
        return getfromB(input,sym)
    end
end

function evalTape(num::Number,tape::SSATape)
    number=tape.number
    if(!hasB(number,num))
        newNode=NumberNode(length(number))
        add(number,newNode=>num)
    end
    return getfromB(number,num)
end

function evalTape(ssaform::SSAForm,tape::SSATape)
    compute=tape.compute
    if(!hasB(compute,ssaform))
        newNode=ComputeNode(length(compute))
        add(compute,newNode=>ssaform)
    end
    return getfromB(compute,ssaform)
end


function evalTape(expr::Expr,tape::SSATape)
    if(expr.head==:call)
        # we need to treat power with special care, right now, we just treat the remaining case
        op=expr.args[1]
        args=[evalTape(arg_,tape) for arg_ in expr.args[2:end]]
        evalTape(SSAForm(op,args),tape)
    else
        error("unexpected Expr $(expr)")
    end    
end

# finally, we need to assign a compute Node to output and with a name
"""
node is assumed to be ComputeNode, but could also be InputNode
"""
function evalTape(node::Node, sym::Symbol,tape::SSATape)
    output_name=tape.output_name
    outputNode=OutputNode(length(output_name))
    tape.output[outputNode]=node
    add(output_name,outputNode=>sym)
end
