# FFSSolver: (F)orward (S)ingle (S)hooting (Solver)
# Solves indirect trajectory optimization problem using a forward
# single shooting based approach.

# NOTE: Likely only supports spacecraft trajectory optimization problems 
# involving a single spacecraft using a 6 element state representation
# with mass (7 total elements). 

struct FSSSolver{NLEF} <: ShootingSolver 

    # Initial guess vector
    initGuess::Vector{Float64}

    # Initial and final boundary conditions
    ics::Vector{Float64}
    fcs::Vector{Float64}

    # Non-linear equations function
    nlEqs::NLEF

    # Data preallocation for function input vectors
    ivNoSTM::Vector{Float64} 
    ivWSTM::Vector{Float64}

end

function FSSSolver(bvpFunc, ICS, FCS)
    # When setting iv(s) need to set indecies coresponding to initGuess!
end

fssNLEqFunc!(F::Union{Nothing,AbstractArray}, J::Union{Nothing,AbstractArray},
             λ::AbstractArray, ivNoSTM, ivWSTM, bvpFuncNoSTM, bvpFuncWSTM, ics, fcs)

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
        ivNoSTM[m+1:n] .= λ

        # Evaluate boundary value function
        zf = bvpFuncNoSTM(ivNoSTM)

        # Fill F (This is not generalized!)
        @views F[1:6] .= zf[1:6] .- fcs[1:6]
        F[7] = zf[14] - fcs[7]
    elseif J !== nothing
        # Initialize full initial condition vector w/ STM
        ivWSTM[m+1:n] .= λ

        # !!! THIS BLOCK OF CODE MAY NOT BE NECESSARY
        ivWSTM[n+1:end] .= 0.0
        for i in 1:n
            ivWSTM[n + (i-1)*n + i] = 1.0
        end
        # !!! END BLOCK

        # Evaluate boundary value function
        zf = bvpFuncWSTM(ivWSTM)

        # Fill F (This is not generalized!)
        if F !== nothing 
            @views F[1:6] .= zf[1:6] .- fcs[1:6]
            F[7] = zf[14] - fcs[7]
        end

        # Fill J
        @views STM = reshape(zf[n+1:end], (n,n))
        @views J .= STM[[1,2,3,4,5,6,14], 8:14] # Not generalized!
    end

    return nothing
end