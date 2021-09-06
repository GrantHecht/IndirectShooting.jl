module IndirectShooting

using NLsolve

# Includes
include("ShootingDataManager.jl")
include("ShootingSolver.jl")
include("FSSSolver.jl")

# Exports
export FSSSolver
export solve!

export initializeData!
export GetInitialGuessConverged
export GetHomotopyConverged
export GetSolution 
export GetHomotopySolutionVector
export GetHomotopyParams
export GetHomotopyConvergenceFlags

end
