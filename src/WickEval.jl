# file to compute the expectation values of a given operator under a non-interacting density matrix parametrized by its single paritcle density matrix (for DA, it becomes the discrete Green function)

# this implement the Wick theorem
# updated version, removing junk code
# this stores a given calculation point
# We should reconstruct the Wick theorem part code.

# We have implement SSATape
# to start, one should create a uniopMap that store
# We add a more flexible way to invoke gc
export evalWick,initUniOpMap
"""
this for SymEngine.jl
"""
function initUniOpMap(input::Dict,ssatape::SSATape)
    uniopMap=Dict{UniOp,Basic}()
    for (k,v) in input
        # print("$(k)->$(v)\n")
        result=evalTape(v,ssatape)
        uniopMap[k]=convert(Basic,result,ssatape)
    end
    uniopMap
end
"""
this is for MathExpr.jl
The idea is to store the result as Node, and avoid parse back and forth the symbolic library and direclty bulid ssatape inplace.
"""
function initUniOpMap(input::Dict{Any,SymExpr},ssatape::SSATape)
    uniopMap=Dict{UniOp,Node}()
    for (k,v) in input
        # print("$(k)->$(v)\n")
        uniopMap[k]=evalTape(v,ssatape)
    end
    uniopMap
end

"""
We should also have a function that evaluate the coefficient only
This may help to avoid unnecessary computation
"""
function evalWick(mop::MultiOp{T,F},ssatape::SSATape) where {T,F}
    s=mop.M
    result=MultiOp(Dict{T,F}(),s)
    for (k,c) in mop.val
        node=ssatape(c)
        value=convert(Basic,node,ssatape)
        result.val[k]=value
    end
    result
end

"""
We should also have a function that evaluate the coefficient only
This may help to avoid unnecessary computation
T=UInt16
k,c=collect(mop.val)[1]
"""
function evalWick(mop::MultiOp{T,SymExpr},ssatape::SSATape) where T
    s=mop.M
    result=MultiOp(Dict{T,SymExpr}(),s)
    engine=__default__engine__map__["default"]
    for (k,c) in mop.val
        node=ssatape(c)
        value=engine(Symbol("$(node)"))
        result.val[k]=value
    end
    result
end


"""
num=160000
# we should group the uop according to the particle number, and evaluate from small to large.
# The number of steps to perform gc may be not much
We should put the num_step_to_gc as an argument
this returns the Node result
"""
function evalWick(mop::MultiOp,uniopMap::Dict{UniOp,Basic},ssatape::SSATape;num_step_to_gc=5000)
    value=Basic("0")
    idx=1
    num_uops=length(mop.val)
    print_idx=1
    num_step_to_print=round(Int,num_uops/1000)
    gc_idx=1
    # 10 is too small
    # num_step_to_gc=round(Int,num_uops/100)
    print("start eval mop with size $(num_uops)\n")
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
    for op_num  in sort(collect(keys(num_to_uop)))
        uops=num_to_uop[op_num]
        for uop in uops
            c=mop.val[uop.val]
            value+=c*evalWick(uop,uniopMap,ssatape)
            if(print_idx>num_step_to_print)
                print("$(idx)/$(num_uops): $(op_num) ops\n")
                print_idx=0
            end
            if(gc_idx>num_step_to_gc)
                print("call gc()\n")
                customized_gc()
                gc_idx=0
            end            
            idx+=1
            print_idx+=1
            gc_idx+=1
        end
    end    
    evalTape(value,ssatape)
end


"""
num=160000
# we should group the uop according to the particle number, and evaluate from small to large.
# The number of steps to perform gc may be not much
We should put the num_step_to_gc as an argument
this returns the Node result
We update to use MathExpr.jl
"""
function evalWick(mop::MultiOp,uniopMap::Dict{UniOp,Node},termMap::Dict{Term,Node}, ssatape::SSATape;num_step_to_gc=5000)
    # value=Basic("0")
    # the result, now, we store all the terms in array
    terms=Array{Node,1}()
    idx=1
    num_uops=length(mop.val)
    print_idx=1
    num_step_to_print=round(Int,num_uops/1000)
    gc_idx=1
    # 10 is too small
    # num_step_to_gc=round(Int,num_uops/100)
    print("start eval mop with size $(num_uops)\n")
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
    # op_num=2
    for op_num  in sort(collect(keys(num_to_uop)))
        uops=num_to_uop[op_num]
        # uop=collect(uops)[1]
        for uop in uops
            c=mop.val[uop.val]
            c_value=ssatape(c,termMap)
            uop_value=evalWick(uop,uniopMap,ssatape)
            if(c_value!=ssatape(0.0) && uop_value!=ssatape(0.0))
                push!(terms,ssatape(SSAForm(:(*),[c_value,uop_value])))
            end            
            if(print_idx>num_step_to_print)
                print("$(idx)/$(num_uops): $(op_num) ops\n")
                print_idx=0
            end
            if(gc_idx>num_step_to_gc)
                print("call gc()\n")
                customized_gc()
                gc_idx=0
            end            
            idx+=1
            print_idx+=1
            gc_idx+=1
        end
    end    
    n_terms=length(terms)
    if(n_terms==0)
        return ssatape(0.0)
    elseif(n_terms==1)
        return terms[1]
    else
        return ssatape(SSAForm(:(+),terms))
    end    
end



"""
evaluate of a uniop, this returns the symbolic result
"""
function evalWick(uop::UniOp, uniopMap::Dict{UniOp,Basic},ssatape::SSATape)
    pos=getOp(uop)              # position of the operators
    num=length(pos)             # number of operators
    zero_=Basic("0")
    one_=Basic("1")
    if(num%2==1)                # we assume the particle conservation
        return zero_
    elseif(num==0)              #  for identity, we normalize the trace of density matrix to 1
        return one_
    elseif(num==2)              # we assume that if a single particle is not in uniopMap, then it is zero, and we don't need to save it
        if(haskey(uniopMap,uop))
            return uniopMap[uop]
        else
            return zero_
        end        
    else
        # prevoius code has too many places to check zero, I think with the current design, these checkings are not necessary.
        # i=2
        if(!haskey(uniopMap,uop))
            value=zero_
            for i in 2:num
                uop_pair=createOp(pos[[1,i]],uop.M)
                uop_pair_value=evalWick(uop_pair,uniopMap,ssatape)
                if(uop_pair_value!=0)
                    pos_rest=copy(pos)
                    deleteat!(pos_rest,[1,i])
                    uop_remain=createOp(pos_rest,uop.M)
                    sign=Basic((-1)^i)
                    uop_remain_value=evalWick(uop_remain,uniopMap,ssatape)
                    value+=uop_pair_value*uop_remain_value*sign
                end
            end
            # now, we evalute expression in ssatape
            value_node=evalTape(value,ssatape)
            value_sym=convert(Basic,value_node,ssatape)
            uniopMap[uop]=value_sym
        end        
        return uniopMap[uop]
    end    
end

"""
evaluate of a uniop, this returns the symbolic result XXXXXX
using MathExpr.jl
now, we return the Node 
"""
function evalWick(uop::UniOp, uniopMap::Dict{UniOp,Node},ssatape::SSATape)
    pos=getOp(uop)              # position of the operators
    num=length(pos)             # number of operators
    # zero_=Basic("0")
    # we just return the Node now
    # one_=Basic("1")
    if(num%2==1)                # we assume the particle conservation
        return ssatape(0.0)
    elseif(num==0)              #  for identity, we normalize the trace of density matrix to 1
        return ssatape(1.0)
    elseif(num==2)              # we assume that if a single particle is not in uniopMap, then it is zero, and we don't need to save it
        if(haskey(uniopMap,uop))
            return uniopMap[uop]
        else
            return ssatape(0.0)
        end        
    else
        # prevoius code has too many places to check zero, I think with the current design, these checkings are not necessary.
        # i=2
        if(!haskey(uniopMap,uop))
            terms=Array{Node,1}()
            # i=2
            for i in 2:num
                uop_pair=createOp(pos[[1,i]],uop.M)
                uop_pair_value=evalWick(uop_pair,uniopMap,ssatape)
                if(uop_pair_value!=ssatape(0.0))
                    pos_rest=copy(pos)
                    deleteat!(pos_rest,[1,i])
                    uop_remain=createOp(pos_rest,uop.M)
                    sign=(-1)^i
                    uop_remain_value=evalWick(uop_remain,uniopMap,ssatape)
                    if(uop_remain_value!=ssatape(0.0))
                        if(sign==-1)
                            # we have update the api, TODO, update the SSAForm call
                            push!(terms,ssatape(SSAForm(:(*),[ssatape(-1),uop_pair_value,uop_remain_value])))
                        else
                            push!(terms,ssatape(SSAForm(:(*),[uop_pair_value,uop_remain_value])))
                        end
                    end                    
                end
            end
            # now, we evalute expression in ssatape
            n_terms=length(terms)
            if(n_terms==0)
                uniopMap[uop]=ssatape(0.0)
            elseif(n_terms==1)
                uniopMap[uop]=terms[1]
            else
                uniopMap[uop]=ssatape(SSAForm(:(+),terms))
            end            
        end        
        return uniopMap[uop]
    end    
end





# old code:



# struct Item{F}
#     name::F
#     value::F
#     order::Int
# end

# """
# for intermediate uniop
# """
# wickReserved="M"

# function wick(mop::MultiOp{T,F},input::Dict{UniOp{T},F},history::Dict{Union{UniOp{T},String},Item{Basic}}) where  {T,F}
#     result=convert(F,0)
#     for (k,c) in mop.val
#         result+=c*wick(UniOp(k,mop.M),input,history)
#     end
#     result
# end

# """
# there are redundant zero nodes in history, try to figure out why.
# """
# function wick(uop::UniOp{T},input::Dict{UniOp{T},F},history::Dict{Union{UniOp{T},String},Item{Basic}}) where  {T,F}
#     numOfOps=getNum(uop)
#     if(numOfOps%2==1)
#         return convert(Basic,0)
#     elseif(numOfOps==0)
#         return convert(Basic,1)
#     elseif(numOfOps==2)
#         if(haskey(input,uop))
#             return input[uop]
#         else
#             return convert(Basic,0)
#         end        
#     else
#         # we need to check whether it is zero or not
#         if(haskey(history,uop))
#             if(history[uop].value!=0)
#                 return history[uop].name
#             else
#                 return convert(Basic,0)
#             end
#         else
#             value=convert(Basic,0)
#             pos=getOp(uop)
#             for i in 2:length(pos)
#                 pos_rest=copy(pos)
#                 deleteat!(pos_rest,[1,i])
#                 factor=wick(createOp([pos[1],pos[i]],uop.M),input,history)
#                 if(factor!=0)
#                     remains=wick(createOp(pos_rest,uop.M),input,history)
#                     if(remains!=0)
#                         if(!haskey(history,remains)||(history[remains].value!=0))
#                             value+=factor*remains*(-1)^i
#                         end                            
#                     end
#                 end
#             end
#             name= convert(Basic,"$(wickReserved)$(uop.val)")
#             history[uop]=Item(name,value,length(history))
#             if(value!=0)
#                 return name
#             else
#                 return convert(Basic,0)
#             end            
#         end        
#     end        
# end

