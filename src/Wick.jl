module Wick
using SSA                       # we move  SSAxxx to a new package
using MathExpr
# include("SSATape.jl")
# include("SSAInterpreter.jl")
# include("SSACompiler.jl")
# include("SSATapeDiff.jl")
include("Memory.jl")
include("Op.jl")
include("WickEval.jl")
include("OpDiff.jl")

end # module
