# FMSSolver: (Forward) (M)ultiple (S)hooting (S)olver

# Abstraction layer to allow for dispatching based on shooting data manager
abstract type AbstractFMSSolver{SDMT <: AbstractShootingDataManager} end

struct FMSSolver{BVPFuncT,SDMT,NLEF} <: AbstractFMSSolver{SDMT}
    # Time span
    tspan::Tuple{Float64, Float64}

    # Number of segments
    nSeg::Int

    # Boundary value function with no STM computation
    bvpFunc::BVPFuncT

    # Shooting data manager
    sdm::SDMT

    # Nonlinear equations 
    nlEqs::NLEF
end

function FMSSolver(tspan, ics, fcs, bvpFuncNoSTM, bvpFuncWSTM;
    nSeg = 4, homotopy = false, homotopyParamVec = [1.0])

    # Get size of state/co-state vector
    m = 7
    n = 2*m

    # Set up ShootingDataManager
    if homotopy == false
        sdm = ShootingDataManager()
    else
        sdm = ShootingHomotopyDataManager()
        SetHomotopyParams!(sdm, homotopyParamVec)
    end
    SetInitGuessData!(sdm, zeros(m + n*(nSeg - 1)))
    SetInitConditions!(sdm, ics)
    SetFinConditions!(sdm, fcs)
    SetPreallocatedVecs!(sdm, n + n*n) 
    Initialize!(sdm)

    # Configure input vector
    inVec = GetInputVec(sdm)
    inVec[1:m]      .= ics 
    inVec[m+1:end]  .= 0.0
    for i in 1:n 
        inVec[n + i + (i-1)*n] = 1.0
    end

    # Configure NLsolve input function
    if homotopy == false
        nlEqs = (F,J,x) -> fmsNLEqFunc!(F, J, x, bvpFuncNoSTM, bvpFuncWSTM,
            tspan, nSeg, sdm)
    else
        nlEqs = (F,J,x,ϵ) -> fmsHomotopyNLEqFunc!(F, J, x, ϵ, bvpFuncNoSTM, bvpFuncWSTM,
            tspan, nSeg, sdm)
    end

    # Construct solver
    FMSSolver(tspan, nSeg, bvpFuncNoSTM, sdm, nlEqs)
end

function initializeData!(solver::FMSSolver, initGuess::AbstractVector)
    m = length(initGuess)
    n = 2*m

    # Compute times
    ts      = zeros(solver.nSeg + 1)
    ts[1]   = solver.tspan[1]
    ts[end] = solver.tspan[2]
    for i in 1:solver.nSeg - 1
        ts[i + 1] = ts[i] + (ts[end] - ts[1]) / solver.nSeg
    end

    # Compute full initial guess vector
    xs      = zeros(m + (solver.nSeg - 1)*n)
    xs[1:m] .= initGuess
    y       = zeros(n)
    for i in 1:solver.nSeg - 1
        if i == 1
            # Construct initial state vector
            y[1:m]  .= GetInitConditions(solver.sdm)
            y[m+1:n].= initGuess

            # Construct time span
            segTspan = (ts[i], ts[i + 1])

            # Perform integration and set next state vector
            idx0 = m + 1
            idxf = idx0 + n - 1
            if solver isa AbstractFMSSolver{ShootingDataManager}
                xs[idx0:idxf] .= solver.bvpFunc(y, segTspan)
            else
                ϵ = GetHomotopyParams(solver.sdm)[1]
                xs[idx0:idxf] .= solver.bvpFunc(y, segTspan, ϵ)
            end
        else
            # Construct initial state vector
            idx0 = m + (i - 2)*n + 1
            idxf = idx0 + n - 1
            y   .= xs[idx0:idxf]

            # Construct time span
            segTspan = (ts[i], ts[i + 1])

            # Perform integration and set next state vector
            idx0 = m + (i - 1)*n + 1
            idxf = idx0 + n - 1
            if solver isa AbstractFMSSolver{ShootingDataManager}
                xs[idx0:idxf] .= solver.bvpFunc(y, segTspan)
            else
                ϵ = GetHomotopyParams(solver.sdm)[1]
                xs[idx0:idxf] .= solver.bvpFunc(y, segTspan, ϵ)
            end
        end
    end
    SetInitGuessData!(solver.sdm, xs)
    SetPreallocatedVecs!(solver.sdm, n + n*n)
    Initialize!(solver.sdm)

    # Set input vector identity matrix portion
    inVec   = GetInputVec(solver.sdm)
    inVec  .= 0.0
    for i in 1:n
        inVec[n + i + (i - 1)*n] = 1.0
    end
    return nothing
end

function solve!(solver::AbstractFMSSolver{ShootingDataManager}; 
    factor = 1.0, ftol = 1e-8, showTrace = true, convergenceAttempts = 4)
    sdm = solver.sdm 
    sol = nlsolve(only_fj!(solver.nlEqs), GetInitGuessData(sdm);
        show_trace = showTrace, factor = factor, ftol = ftol)
    solver.sol .= sol.zero
    sdm.cflag = (sol.residual_norm < ftol*100.0)
    return nothing
end

function solve!(solver::AbstractFMSSolver{ShootingHomotopyDataManager};
    factor = 1.0, ftol = 1e-8, showTrace = true, convergenceAttepts = 4)
    sdm         = solver.sdm
    ϵs          = GetHomotopyParams(sdm)
    sols        = GetHomotopySolutionVector(sdm)
    cflags      = GetHomotopyConvergenceFlags(sdm)
    lastConvIdx = 1
    cFlag       = false
    @inbounds for i in eachindex(ϵs)
        # Get homotopy parameter and solve
        ϵ = ϵs[i]
        if i == 1
            sol = nlsolve(only_fj!((F,J,x)->solver.nlEqs(F,J,x,ϵ)),
                GetInitGuessData(sdm); show_trace = showTrace, 
                factor = factor, ftol = ftol)
        else
            println("Shooting for ϵ = " * string(ϵ))
            if cFlag == true
                lastConvIdx = i - 1
            end
            sol = nlsolve(only_fj!((F,J,x)->solver.nlEqs(F,J,x,ϵ)), sols[lastConvIdx];
                show_trace = showTrace, factor = factor, ftol = ftol) 
        end

        # Save solution
        cFlag       = (sol.residual_norm < ftol*100.0)
        cflags[i]   = cFlag
        sols[i]    .= sol.zero

        # Edge case catch: If final homotopy step fails, attempt to correct
        # be reducing steps.
        if i == length(ϵs) && cflags[end] == false && cflags[end-1] == true
            println("Last step failed. Attempting to correct...")
            attempts  = 0
            success   = false
            localSol  = Vector{Float64}(undef, length(sols[1]))
            localSol .= sols[end - 1]
            localSolϵ = ϵs[end - 1]
            while success == false && attempts <= convergenceAttepts
                attempts    += 1
                numSteps    = 2 + attempts
                for ϵ in LinRange(localSolϵ, ϵs[end], numSteps)[2:end]
                    println("Shooting for ϵ = " * string(ϵ))
                    sol = nlsolve(only_fj!((F,J,x)->solver.nlEqs(F,J,x,ϵ)),
                        localSol; show_trace = showTrace, factor = factor, ftol = ftol)
                    cFlag = (sol.residual_norm < ftol*100.0)
                    if cFlag == false
                        break
                    else
                        localSol .= sol.zero
                        localSolϵ = ϵ
                    end
                end
                if cFlag == true
                    success = true
                    sols[end] .= localSol
                    cflags[end] = true
                end
            end
        end

        # Check if convergence has stalled
        if cFlag == false && (i - lastConvIdx) >= convergenceAttepts
            println("Convergence failed during continuation.")
            break
        end
    end
end

# Nonlinear functions
function fmsNLEqFunc!(F::Union{Nothing,AbstractArray}, J::Union{Nothing,AbstractArray},
    x::AbstractArray, bvpFuncNoSTM, bvpFuncWSTM, tspan, nSeg,
    sdm::ShootingDataManager)

    # Get data from shooting data manager
    z0 = GetInputVec(sdm)
    zf = GetOutputVec(sdm)
    ic = GetInitConditions(sdm) 
    fc = GetFinConditions(sdm)

    # Compute times
    ts      = zeros(nSeg + 1)
    ts[1]   = tspan[0]
    ts[end] = tspan[1]
    for i in 1:nSeg - 1
        ts[i + 1] = ts[i] + (ts[end] - ts[1]) / nSeg
    end

    # Fill J with zeros if necesary
    if !(J === nothing)
        fill!(J, 0.0)
    end

    # Set helper vars
    m = 7
    n = 14

    # Evaluate functions
    if J === nothing && F !== nothing
        for i in 1:nSeg
            # Fill initial state vector
            if i == 1
                @views z0[1:m]  .= ic
                @views z0[m+1:n].= x[1:m]
            else
                idx0 = m + (i - 2)*n + 1
                idxf = m + (i - 1)*n
                @views z0[1:n]  .= x[idx0:idxf]
            end

            # Compute segment tspan
            segTspan = (ts[i], ts[i + 1])

            # Evaluate function
            @views zf[1:n] .= bvpFuncNoSTM(z0[1:n], segTspan)

            # Fill constraint vector
            if i == nSeg
                idx0 = (i - 1)*n + 1
                idxf = idx0 + 5
                @views F[idx0:idxf] .= zf[1:6] .- fc[1:6]
                F[end] = zf[14] - fc[7]
            else
                idx0 = (i - 1)*n + 1
                idxf = i*n
                @views F[idx0:idxf] .= zf[1:n] .- x[m + idx0:m + idxf] 
            end
        end
    elseif J!== nothing
        for i in 1:nSeg
            # Fill initial state vector
            if i == 1
                @views z0[1:m]  .= ic
                @views z0[m+1:n].= x[1:m]
            else
                idx0 = m + (i - 2)*n + 1
                idxf = m + (i - 1)*n
                @views z0[1:n]  .= x[idx0:idxf]
            end

            # Compute segment tspan
            segTspan = (ts[i], ts[i + 1])

            # Evaluate function
            zf .= bvpFuncWSTM(z0, segTspan)

            # Fill constraint vector and Jacobian matrix
            @views STM = reshape(zf[n+1:end], (n,n))
            if i == nSeg
                xidx0 = m + (i - 2)*n + 1
                xidxf = m + (i - 1)*n
                fidx0 = (i - 1)*n + 1
                fidxf = idx0 + 5

                if F !== nothing
                    @views F[fidx0:fidxf] .= zf[1:6] .- fc[1:6]
                    F[end] = zf[14] - fc[7]
                end
                @views J[idx0:idxf+1,xidx0:xidxf] .= STM[[1,2,3,4,5,6,14], :]
            else
                fidx0 = (i - 1)*n + 1
                fidxf = i*n
            
                if F !== nothing
                    @views F[fidx0:fidxf] .= zf[1:n] .- x[m + idx0:m + idxf] 
                end

                if i == 1
                    @views J[1:14,1:7] .= STM[:, [8:14]]
                else
                    xidx0 = m + (i - 2)*n + 1
                    xidxf = m + (i - 1)*n
                    @views J[fidx0:fidxf,xidx0:xidxf] .= STM
                end
            end
        end
    end
    return nothing
end

function fmsHomotopyNLEqFunc!(F::Union{Nothing,AbstractArray}, J::Union{Nothing,AbstractArray},
    x::AbstractArray, ϵ, bvpFuncNoSTM, bvpFuncWSTM, tspan, nSeg,
    sdm::ShootingHomotopyDataManager)

    # Get data from shooting data manager
    z0 = GetInputVec(sdm)
    zf = GetOutputVec(sdm)
    ic = GetInitConditions(sdm) 
    fc = GetFinConditions(sdm)

    # Compute times
    ts      = zeros(nSeg + 1)
    ts[1]   = tspan[1]
    ts[end] = tspan[2]
    for i in 1:nSeg - 1
        ts[i + 1] = ts[i] + (ts[end] - ts[1]) / nSeg
    end

    # Fill J with zeros if necesary
    if !(J === nothing)
        fill!(J, 0.0)
    end

    # Set helper vars
    m = 7
    n = 14

    # Evaluate functions
    if J === nothing && F !== nothing
        for i in 1:nSeg
            # Fill initial state vector
            if i == 1
                @views z0[1:m]  .= ic
                @views z0[m+1:n].= x[1:m]
            else
                idx0 = m + (i - 2)*n + 1
                idxf = m + (i - 1)*n
                @views z0[1:n]  .= x[idx0:idxf]
            end

            # Compute segment tspan
            segTspan = (ts[i], ts[i + 1])

            # Evaluate function
            @views zf[1:n] .= bvpFuncNoSTM(z0[1:n], segTspan, ϵ)

            # Fill constraint vector
            if i == nSeg
                idx0 = (i - 1)*n + 1
                idxf = idx0 + 5
                @views F[idx0:idxf] .= zf[1:6] .- fc[1:6]
                F[end] = zf[14] - fc[7]
            else
                idx0 = (i - 1)*n + 1
                idxf = i*n
                @views F[idx0:idxf] .= zf[1:n] .- x[m + idx0:m + idxf] 
            end
        end
    elseif J!== nothing
        for i in 1:nSeg
            # Fill initial state vector
            if i == 1
                @views z0[1:m]  .= ic
                @views z0[m+1:n].= x[1:m]
            else
                idx0 = m + (i - 2)*n + 1
                idxf = m + (i - 1)*n
                @views z0[1:n]  .= x[idx0:idxf]
            end

            # Compute segment tspan
            segTspan = (ts[i], ts[i + 1])

            # Evaluate function
            zf .= bvpFuncWSTM(z0, segTspan, ϵ)

            # Fill constraint vector and Jacobian matrix
            @views STM = reshape(zf[n+1:end], (n,n))
            if i == nSeg
                xidx0 = m + (i - 2)*n + 1
                xidxf = m + (i - 1)*n
                fidx0 = (i - 1)*n + 1
                fidxf = fidx0 + 5

                if F !== nothing
                    @views F[fidx0:fidxf] .= zf[1:6] .- fc[1:6]
                    F[end] = zf[14] - fc[7]
                end
                @views J[fidx0:fidxf+1,xidx0:xidxf] .= STM[[1,2,3,4,5,6,14], :]
            else
                xidx0 = m + (i - 2)*n + 1
                xidxf = m + (i - 1)*n
                fidx0 = (i - 1)*n + 1
                fidxf = i*n
            
                if F !== nothing
                    @views F[fidx0:fidxf] .= zf[1:n] .- x[xidx0 + n:xidxf + n] 
                end

                # STM terms 
                if i == 1
                    @views J[1:14,1:7] .= STM[:, 8:14]
                else
                    @views J[fidx0:fidxf,xidx0:xidxf] .= STM
                end

                # Identity terms
                xidx0 = m + (i - 1)*n + 1
                xidxf = m + i*n
                for j in 0:6
                    J[fidx0 + j,xidx0 + j] = -1.0
                end
            end
        end
    end
    return nothing
end

function GetInitialGuessConverged(fmss::FMSSolver)
    return GetInitialGuessConverged(fmss.sdm)
end

function GetHomotopyConverged(fmss::AbstractFMSSolver{ShootingHomotopyDataManager})
    return GetHomotopyConverged(fmss.sdm)
end

function GetSolution(fmss::AbstractFMSSolver{ShootingDataManager})
    return GetSolution(fmss.sdm)
end

function GetHomotopySolutionVector(fmss::AbstractFMSSolver{ShootingHomotopyDataManager})
    return GetHomotopySolutionVector(fmss.sdm)
end

function GetHomotopyParams(fmss::AbstractFMSSolver{ShootingHomotopyDataManager})
    return GetHomotopyParams(fmss.sdm)
end

function GetHomotopyConvergenceFlags(fmss::AbstractFMSSolver{ShootingHomotopyDataManager})
    return GetHomotopyConvergenceFlags(fmss.sdm)
end