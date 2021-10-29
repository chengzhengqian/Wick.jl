# file to compute the expectation values of a given operator under a non-interacting density matrix parametrized by its single paritcle density matrix (for DA, it becomes the discrete Green function)

# this implement the Wick theorem
# updated version, removing junk code
# this stores a given calculation point
# We should reconstruct the Wick theorem part code.

# We have implement SSATape
# to start, one should create a uniopMap that store
# We add a more flexible way to invoke gc
export evalWick,initUniOpMap

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

