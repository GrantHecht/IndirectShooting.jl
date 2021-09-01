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

function solve!(solver::AbstractFSSSolver{ShootingDataManager})
    sdm = solver.sdm
    sol = nlsolve(only_fj!(solver.nlEqs), GetInitGuessData(sdm); show_trace = true)
end

function solve!(solver::AbstractFSSSolver{ShootingHomotopyDataManager})
    sdm = solver.sdm
    ϵs = GetHomotopyParams(solver.sdm)
    sols = GetHomotopySolutionVector(solver.sdm)
    @inbounds for i in 1:length(ϵs)
        # Get homotopy parameter and solve
        ϵ = ϵs[i]
        if i == 1
            sol = nlsolve(only_fj!((F,J,λ)->solver.nlEqs(F,J,λ,ϵ)), 
                GetInitGuessData(sdm); show_trace = true)
        else
            sol = nlsolve(only_fj!((F,J,λ)->solver.nlEqs(F,J,λ,ϵ)),
                sols[i-1]; show_trace = true)
        end
        # Save solution
        sols[i] .= sol.zero
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
        @views zf[1:n] .= bvpFuncNoSTM(z0[1:n])

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