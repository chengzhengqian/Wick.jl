export opDiff
"""
result should be a Mop
"""
function opDiff(target::UniOp,arg::UniOp)
    s=target.M
    sign,result=opDiff(getOp(target),getOp(arg))
    MultiOp(Basic(sign),UniOp(result,s))
end

"""
to implement this, we just need the representation of a UniOp as a array

"""
function opDiff(target_::Array,arg_::Int)
    idx=findall(x->x==arg_,target_)
    if(length(idx)>0)
        idx_=idx[1]
        sign=(-1)^(idx_-1)
        target_=copy(target_)
        deleteat!(target_,idx_)
        return (sign,target_)
    else
        return (0,Array{Int,1}())
    end    
end

"""
assume len(args_)>0
"""
function opDiff(target_::Array,args_::Array)
    sign=1
    for arg in args_
        sign_, target_=opDiff(target_,arg)
        if(sign_==0)
            return (0,Array{Int,1}())
        end
        sign*=sign_
    end
    return (sign,target_)
end


function opDiff(target::MultiOp{T,F},arg::UniOp) where {T,F}
    s=target.M
    result=MultiOp(Dict{T,F}(),s)
    for (k,c) in target.val
        addMOpMOp!(result,c*opDiff(UniOp(k,s),arg))
    end
    result
end

function opDiff(target::MultiOp{T,F}, args...) where {T,F}
    for arg_ in args
        target=opDiff(target,arg_)
    end
    target
end
