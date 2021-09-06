# FFSSolver: (F)orward (S)ingle (S)hooting (Solver)
# Solves indirect trajectory optimization problem using a forward
# single shooting based approach.

# NOTE: Likely only supports spacecraft trajectory optimization problems 
# involving a single spacecraft using a 6 element state representation
# with mass (7 total elements). 

# Abstraction layer to allow for disptaching based on ShootingDataManager
abstract type AbstractFSSSolver{SDMT <: AbstractShootingDataManager} end

struct FSSSolver{SDMT,NLEF} <: AbstractFSSSolver{SDMT}

    # Shooting Data Manager 
    sdm::SDMT

    # Non-linear equations function
    nlEqs::NLEF

end

# Constructor
function FSSSolver(initGuess, ics, fcs, bvpFuncNoSTM, bvpFuncWSTM;
    homotopy = false, homotopyParamVec = [1.0])

    # Get size of variables vector
    m = length(initGuess)
    n = 2*m 

    # Set up ShootingDataManager
    if homotopy == false
        sdm = ShootingDataManager()
    else
        sdm = ShootingHomotopyDataManager()
        SetHomotopyParams!(sdm, homotopyParamVec)
    end
    SetInitGuessData!(sdm, initGuess)
    SetInitConditions!(sdm, ics)
    SetFinConditions!(sdm, fcs)
    SetPreallocatedVecs!(sdm, n + n*n)
    Initialize!(sdm)

    # Configuring input vector
    inVec = GetInputVec(sdm)
    inVec[1:m]      .= ics
    inVec[m+1:n]    .= initGuess
    inVec[n+1:end]  .= 0.0
    for i in 1:n
        inVec[n + i + (i-1)*n] = 1.0
    end

    # Configure NLsolve.jl input function
    if homotopy == false
        nlEqs = (F,J,λ) -> fssNLEqFunc!(F, J, λ, bvpFuncNoSTM, bvpFuncWSTM, sdm)
    else 
        nlEqs = (F,J,λ,ϵ) -> fssHomotopyNLEqFunc!(F, J, λ, ϵ, bvpFuncNoSTM, bvpFuncWSTM, sdm)
    end

    FSSSolver(sdm, nlEqs)
end

function FSSSolver(ics, fcs, bvpFuncNoSTM, bvpFuncWSTM;
    homotopy = false, homotopyParamVec = [1.0])

    # Set up ShootingDataManager
    if homotopy == false
        sdm = ShootingDataManager()
    else
        sdm = ShootingHomotopyDataManager()
        SetHomotopyParams!(sdm, homotopyParamVec)
    end
    SetInitConditions!(sdm, ics)
    SetFinConditions!(sdm, fcs)

    # Configure NLsolve.jl input function
    if homotopy == false
        nlEqs = (F,J,λ) -> fssNLEqFunc!(F, J, λ, bvpFuncNoSTM, bvpFuncWSTM, sdm)
    else 
        nlEqs = (F,J,λ,ϵ) -> fssHomotopyNLEqFunc!(F, J, λ, ϵ, bvpFuncNoSTM, bvpFuncWSTM, sdm)
    end

    FSSSolver(sdm, nlEqs)
end


function initializeData!(solver::FSSSolver, initGuess::AbstractVector)
    m = length(initGuess)
    n = 2*m
    SetInitGuessData!(solver.sdm, initGuess)
    SetPreallocatedVecs!(solver.sdm, n + n*n)
    Initialize!(solver.sdm)
    inVec = GetInputVec(solver.sdm)
    inVec[1:m]      .= solver.sdm.ics 
    inVec[m+1:n]    .= initGuess
    inVec[n+1:end]  .= 0.0
    for i in 1:n
        inVec[n + i + (i-1)*n] = 1.0
    end

    return nothing
end

function solve!(solver::AbstractFSSSolver{ShootingDataManager}; factor = 1.0, ftol = 1e-8)
    sdm = solver.sdm
    sol = nlsolve(only_fj!(solver.nlEqs), GetInitGuessData(sdm);
        show_trace = true, factor = factor, ftol =ftol)
    solver.sol .= sol.zero
    sdm.cflag = (sol.residual_norm < ftol*100.0)
    return nothing
end

function solve!(solver::AbstractFSSSolver{ShootingHomotopyDataManager}; factor = 1.0, ftol = 1e-8, convergenceAttempts = 4)
    sdm         = solver.sdm
    ϵs          = GetHomotopyParams(sdm)
    sols        = GetHomotopySolutionVector(sdm)
    cflags      = GetHomotopyConvergenceFlags(sdm)
    lastConvIdx = 1
    cFlag       = false
    @inbounds for i in 1:length(ϵs)
        # Get homotopy parameter and solve
        ϵ = ϵs[i]
        if i == 1
            sol = nlsolve(only_fj!((F,J,λ)->solver.nlEqs(F,J,λ,ϵ)), 
                GetInitGuessData(sdm); show_trace = true, factor = factor, ftol = ftol)
        else
            println("Shooting for ϵ = " * string(ϵ))
            if cFlag == true
                lastConvIdx = i - 1
            end
            sol = nlsolve(only_fj!((F,J,λ)->solver.nlEqs(F,J,λ,ϵ)),
                sols[lastConvIdx]; show_trace = true, factor = factor, ftol = ftol)
        end

        # Save solution
        cFlag = (sol.residual_norm < ftol*100.0)
        cflags[i] = cFlag
        sols[i] .= sol.zero

        # Edge case catch: If final homotopy step fails, attempt to correct
        # by reducing steps.
        if i == length(ϵs) && cflags[end] == false && cflags[end-1] == true
            println("Last step failed. Attempting to correct...")
            attempts    = 0
            success     = false
            localSol    = Vector{Float64}(undef, length(sols[1]))
            while success == false && attempts <= convergenceAttempts
                attempts    += 1
                numSteps    = 2 + attempts
                localSol    .= sols[end-1]
                for ϵ in LinRange(ϵs[end-1], ϵs[end], numSteps)[2:end] # This is allocating twice. Might want to improve but hitting this edge case is rare...
                    println("Shooting for ϵ = " * string(ϵ))
                    sol = nlsolve(only_fj!((F,J,λ)->solver.nlEqs(F,J,λ,ϵ)),
                        localSol; show_trace = true, factor = factor, ftol = ftol)
                    cFlag = (sol.residual_norm < ftol*100.0)
                    if cFlag == false
                        break
                    else
                        localSol .= sol.zero
                    end
                end
                if cFlag == true
                    success = true
                    sols[end] .= localSol
                end
            end
        end

        # Check if convergence has stalled 
        if cFlag == false && (i - lastConvIdx) >= convergenceAttempts
            println("Convergence failed during continuation.")
            break
        end
    end
end

# NLsolve.jl function
function fssNLEqFunc!(F::Union{Nothing,AbstractArray}, J::Union{Nothing,AbstractArray},
                      λ::AbstractArray, bvpFuncNoSTM, bvpFuncWSTM, 
                      sdm::ShootingDataManager)

    # Get data from ShootingDataManager
    z0 = GetInputVec(sdm)
    zf = GetOutputVec(sdm)
    fc = GetFinConditions(sdm)

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
        z0[m+1:n] .= λ

        # Evaluate boundary value function
        @views zf[1:n] .= bvpFuncNoSTM(z0[1:n])

        # Fill F (This is not generalized!)
        @views F[1:6] .= zf[1:6] .- fc[1:6]
        F[7] = zf[14] - fc[7]
    elseif J !== nothing
        # Initialize full initial condition vector w/ STM
        z0[m+1:n] .= λ
n
        # Evaluate boundary value function with STM
        zf .= bvpFuncWSTM(z0)

        # Fill F (This is not generalized!)
        if F !== nothing 
            @views F[1:6] .= zf[1:6] .- fc[1:6]
            F[7] = zf[14] - fc[7]
        end

        # Fill J
        @views STM = reshape(zf[n+1:end], (n,n))
        @views J .= STM[[1,2,3,4,5,6,14], 8:14] # Not generalized!
    end

    return nothing
end

# NLsolve.jl function
function fssHomotopyNLEqFunc!(F::Union{Nothing,AbstractArray}, J::Union{Nothing,AbstractArray},
                              λ::AbstractArray, ϵ, bvpFuncNoSTM, bvpFuncWSTM, 
                              sdm::ShootingHomotopyDataManager)

    # Get data from ShootingDataManager
    z0 = GetInputVec(sdm)
    zf = GetOutputVec(sdm)
    fc = GetFinConditions(sdm)

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
        z0[m+1:n] .= λ

        # Evaluate boundary value function
        @views zf[1:n] .= bvpFuncNoSTM(z0[1:n], ϵ)

        # Fill F (This is not generalized!)
        @views F[1:6] .= zf[1:6] .- fc[1:6]
        F[7] = zf[14] - fc[7]
    elseif J !== nothing
        # Initialize full initial condition vector w/ STM
        z0[m+1:n] .= λ
n
        # Evaluate boundary value function with STM
        zf .= bvpFuncWSTM(z0, ϵ)

        # Fill F (This is not generalized!)
        if F !== nothing 
            @views F[1:6] .= zf[1:6] .- fc[1:6]
            F[7] = zf[14] - fc[7]
        end

        # Fill J
        @views STM = reshape(zf[n+1:end], (n,n))
        @views J .= STM[[1,2,3,4,5,6,14], 8:14] # Not generalized!
    end

    return nothing
end

function GetInitialGuessConverged(fsss::FSSSolver)
    return GetInitialGuessConverged(fsss.sdm)
end

function GetHomotopyConverged(fsss::AbstractFSSSolver{ShootingHomotopyDataManager})
    return GetHomotopyConverged(fsss.sdm)
end

function GetSolution(fsss::AbstractFSSSolver{ShootingDataManager})
    return GetSolution(fsss.sdm)
end

function GetHomotopySolutionVector(fsss::AbstractFSSSolver{ShootingHomotopyDataManager})
    return GetHomotopySolutionVector(fsss.sdm)
end

function GetHomotopyParams(fsss::AbstractFSSSolver{ShootingHomotopyDataManager})
    return GetHomotopyParams(fsss.sdm)
end

function GetHomotopyConvergenceFlags(fsss::AbstractFSSSolver{ShootingHomotopyDataManager})
    return GetHomotopyConvergenceFlags(fsss.sdm)
end