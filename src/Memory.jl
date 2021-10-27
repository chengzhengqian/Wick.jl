# we add a customized gc which will be triggered when the memory increased by a given amount compared with the previous record value
# after some trial. It is simply just call gc in any case.
function get_memory_usage()
    parse(Int,split(String(read(`ps -p $(getpid()) -o rss`)),"\n")[2])
    # out = read(`top -bn1 -p $(getpid())`, String)
    # res=split(split(out,  "\n")[end-1])[6]
    # if(res[end]=='g')
    #     round(Int,parse(Float64,res[1:(end-1)])*1000000)
    # elseif(res[end]=='m')
    #     round(Int,parse(Float64,res[1:(end-1)])*1000)
    # else
    #     parse(Int,res)
    # end    
end

# get_memory_usage()
# in kB
__previous_memory_usage__=0
# ~500MB
__memory_increase_threshold__=100000


function customized_gc()
    global __previous_memory_usage__, __memory_increase_threshold__
    current_memory_usage=get_memory_usage()
    # if(current_memory_usage<__previous_memory_usage__)
    #     print("$(__previous_memory_usage__)->$(current_memory_usage); does not perform gc, update  __previous_memory_usage__\n")
    #     __previous_memory_usage__=current_memory_usage
    # elseif (current_memory_usage<__previous_memory_usage__+__memory_increase_threshold__)
    #     print("$(__previous_memory_usage__)->$(current_memory_usage); does not perform gc, keep  __previous_memory_usage__\n")
    # else
    print("$(__previous_memory_usage__)->$(current_memory_usage);  perform gc\n")
    [GC.gc() for _ in 1:2]
    __previous_memory_usage__=get_memory_usage()
    print("after gc: $(current_memory_usage)->$(__previous_memory_usage__)\n")
    # end    
end

