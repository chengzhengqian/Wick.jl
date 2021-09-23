# file to compute the expectation values of a given operator under a non-interacting density matrix parametrized by its single paritcle density matrix (for DA, it becomes the discrete Green function)

# this implement the Wick theorem
# updated version, removing junk code
# this stores a given calculation point
# We should reconstruct the Wick theorem part code.




struct Item{F}
    name::F
    value::F
    order::Int
end

"""
for intermediate uniop
"""
wickReserved="M"

function wick(mop::MultiOp{T,F},input::Dict{UniOp{T},F},history::Dict{Union{UniOp{T},String},Item{Basic}}) where  {T,F}
    result=convert(F,0)
    for (k,c) in mop.val
        result+=c*wick(UniOp(k,mop.M),input,history)
    end
    result
end

"""
there are redundant zero nodes in history, try to figure out why.
"""
function wick(uop::UniOp{T},input::Dict{UniOp{T},F},history::Dict{Union{UniOp{T},String},Item{Basic}}) where  {T,F}
    numOfOps=getNum(uop)
    if(numOfOps%2==1)
        return convert(Basic,0)
    elseif(numOfOps==0)
        return convert(Basic,1)
    elseif(numOfOps==2)
        if(haskey(input,uop))
            return input[uop]
        else
            return convert(Basic,0)
        end        
    else
        # we need to check whether it is zero or not
        if(haskey(history,uop))
            if(history[uop].value!=0)
                return history[uop].name
            else
                return convert(Basic,0)
            end
        else
            value=convert(Basic,0)
            pos=getOp(uop)
            for i in 2:length(pos)
                pos_rest=copy(pos)
                deleteat!(pos_rest,[1,i])
                factor=wick(createOp([pos[1],pos[i]],uop.M),input,history)
                if(factor!=0)
                    remains=wick(createOp(pos_rest,uop.M),input,history)
                    if(remains!=0)
                        if(!haskey(history,remains)||(history[remains].value!=0))
                            value+=factor*remains*(-1)^i
                        end                            
                    end
                end
            end
            name= convert(Basic,"$(wickReserved)$(uop.val)")
            history[uop]=Item(name,value,length(history))
            if(value!=0)
                return name
            else
                return convert(Basic,0)
            end            
        end        
    end        
end

