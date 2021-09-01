module IndirectShooting

using NLsolve

# Includes
include("ShootingDataManager.jl")
include("ShootingSolver.jl")
include("FSSSolver.jl")

# Exports
export FSSSolver
export solve!

end
