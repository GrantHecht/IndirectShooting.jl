abstract type AbstractShootingDataManager end

# ShootingDataManager: Standard shooting data manager for solving an 
# indirect trajectory optimization problem. Does not include any data
# for continuation
struct ShootingDataManager <: AbstractShootingDataManager
    # Vector containing initial guess data
    initGuessData::Vector{Float64}

    # Initial and final boundary conditions 
    ics::Vector{Float64}
    fcs::Vector{Float64}

    # Preallocated function input vectors
    inVec::Vector{Float64}
    outVec::Vector{Float64}
end

function ShootingDataManager() 
    ShootingDataManager(Vector{Float64}(undef, 0), 
                        Vector{Float64}(undef, 0),
                        Vector{Float64}(undef, 0),
                        Vector{Float64}(undef, 0),
                        Vector{Float64}(undef, 0))
end

# ShootingHomotopyDataManager: Shooting data manager when performing
# homotopic continuation
struct ShootingHomotopyDataManager <: AbstractShootingDataManager
    # Vector containing initial guess data
    initGuessData::Vector{Float64}

    # Initial and final boundary conditions 
    ics::Vector{Float64}
    fcs::Vector{Float64}

    # Preallocated function input vectors
    inVec::Vector{Float64}
    outVec::Vector{Float64}

    # Continuation parameter vector
    ϵs::Vector{Float64}

    # continuation solution vector
    sols::Vector{Vector{Float64}}
end

function ShootingHomotopyDataManager() 
    ShootingHomotopyDataManager(Vector{Float64}(undef, 0), 
                        Vector{Float64}(undef, 0),
                        Vector{Float64}(undef, 0),
                        Vector{Float64}(undef, 0),
                        Vector{Float64}(undef, 0),
                        Vector{Float64}(undef, 0),
                        Vector{Vector{Float64}}(undef, 0))
end

function SetInitGuessData!(sdm::AbstractShootingDataManager, igv::AbstractVector)
    resize!(sdm.initGuessData, length(igv))
    sdm.initGuessData .= igv
    return nothing
end

function SetInitConditions!(sdm::AbstractShootingDataManager, icv::AbstractVector)
    resize!(sdm.ics, length(icv))
    sdm.ics .= icv 
    return nothing 
end

function SetFinConditions!(sdm::AbstractShootingDataManager, fcv::AbstractVector)
    resize!(sdm.fcs, length(fcv))
    sdm.fcs .= fcv 
    return nothing
end

function SetPreallocatedVecs!(sdm::AbstractShootingDataManager, size::Int)
    resize!(sdm.inVec, size)
    resize!(sdm.outVec, size)
    return nothing
end

function SetPreallocatedVecs!(sdm::AbstractShootingDataManager, inVec::AbstractVector)
    n = length(inVec)
    resize!(sdm.inVec, n)
    resize!(sdm.outVec, n)
    sdm.inVec .= inVec 
    return nothing 
end

function GetInitGuessData(sdm::AbstractShootingDataManager)
    return sdm.initGuessData
end

function GetInitConditions(sdm::AbstractShootingDataManager)
    return sdm.ics
end

function GetFinConditions(sdm::AbstractShootingDataManager)
    return sdm.fcs
end

function GetInputVec(sdm::AbstractShootingDataManager)
    return sdm.inVec
end

function GetOutputVec(sdm::AbstractShootingDataManager)
    return sdm.outVec 
end

function SetHomotopyParams!(sdm::ShootingHomotopyDataManager, ϵs::AbstractVector)
    n = length(ϵs)
    resize!(sdm.ϵs, n)
    sdm.ϵs .= ϵs
    return nothing
end

function GetHomotopyParams(sdm::ShootingHomotopyDataManager)
    return sdm.ϵs
end

function Initialize!(sdm::AbstractShootingDataManager)
    return nothing 
end

function Initialize!(sdm::ShootingHomotopyDataManager)
    # Need to size sols vector of vectors
    n = length(sdm.ϵs)
    m = length(sdm.initGuessData)
    for i in 1:n
        push!(sdm.sols, Vector{Float64}(undef, m))
    end
    return nothing 
end

function GetHomotopySolutionVector(sdm::ShootingHomotopyDataManager)
    return sdm.sols
end
