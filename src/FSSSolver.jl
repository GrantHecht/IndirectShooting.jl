# FFSSolver: (F)orward (S)ingle (S)hooting (Solver)
# Solves indirect trajectory optimization problem using a forward
# single shooting based approach.

# NOTE: Likely only supports spacecraft trajectory optimization problems 
# involving a single spacecraft using a 6 element state representation
# with mass (7 total elements). 

struct FSSSolver{NLEF} <: ShootingSolver 

    # Shooting Data Manager 
    sdm::ShootingDataManager

    # Non-linear equations function
    nlEqs::NLEF

end

# Constructor
function FSSSolver(initGuess, ics, fcs, bvpFuncNoSTM, bvpFuncWSTM)

    # Get size of variables vector
    m = length(initGuess)
    n = 2*m 

    # Set up ShootingDataManager
    sdm = ShootingDataManager()
    SetInitGuessData!(sdm, initGuess)
    SetInitConditions!(sdm, ics)
    SetFinConditions!(sdm, fcs)
    SetPreallocatedVecs!(sdm, n + n*n)

    # Configuring input vector
    sdm.inVec[1:m]      .= ics
    sdm.inVec[m+1:n]    .= initGuess
    sdm.inVec[n+1:end]  .= 0.0
    for i in 1:n
        sdm.inVec[n + i + (i-1)*n] = 1.0
    end

    # Configure NLsolve.jl input function
    nlEqs = only_fj!((F,J,λ) -> fssNLEqFunc!(F, J, λ, bvpFuncNoSTM, bvpFuncWSTM, sdm))

    FSSSolver(sdm, nlEqs)
end

function solve!(solver::FSSSolver)
    sol = nlsolve(solver.nlEqs, solver.sdm.initGuessData; show_trace = true)
end

# NLsolve.jl function
function fssNLEqFunc!(F::Union{Nothing,AbstractArray}, J::Union{Nothing,AbstractArray},
                      λ::AbstractArray, bvpFuncNoSTM, bvpFuncWSTM, 
                      sdm::ShootingDataManager)

    # Fill J with zeros if necessary
    if !(J === nothing)
        fill!(J, 0.0)
    end

    # Added generality to hopfully make current revision applicable to 
    # more problems
    m = length(λ)
    n = 2*m

    # Evaluate boundary value problem and fill data arrays
    if J === nothing && F !== nothing
        # Initialize full initial condition vector w/o STM
        sdm.inVec[m+1:n] .= λ

        # Evaluate boundary value function
        @views sdm.outVec[1:n] .= bvpFuncNoSTM(sdm.inVec[1:n])

        # Fill F (This is not generalized!)
        @views F[1:6] .= sdm.outVec[1:6] .- sdm.fcs[1:6]
        F[7] = sdm.outVec[14] - sdm.fcs[7]
    elseif J !== nothing
        # Initialize full initial condition vector w/ STM
        sdm.inVec[m+1:n] .= λ

        # !!! THIS BLOCK OF CODE MAY NOT BE NECESSARY
        sdm.inVec[n+1:end] .= 0.0
        for i in 1:n
            sdm.inVec[n + (i-1)*n + i] = 1.0
        end
        # !!! END BLOCK

        # Evaluate boundary value function
        sdm.outVec .= bvpFuncWSTM(sdm.inVec)

        # Fill F (This is not generalized!)
        if F !== nothing 
            @views F[1:6] .= sdm.outVec[1:6] .- sdm.fcs[1:6]
            F[7] = sdm.outVec[14] - sdm.fcs[7]
        end

        # Fill J
        @views STM = reshape(sdm.outVec[n+1:end], (n,n))
        @views J .= STM[[1,2,3,4,5,6,14], 8:14] # Not generalized!
    end

    return nothing
end