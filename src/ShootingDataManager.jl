struct ShootingDataManager
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
    undefVec = Vector{Float64}(undef, 0)
    ShootingDataManager(undefVec, undefVec, undefVec, undefVec, undefVec)
end

function SetInitGuessData!(sdm::ShootingDataManager, igv::AbstractVector)
    resize!(sdm.initGuessData, length(igv))
    sdm.initGuessData .= igv
    return nothing
end

function SetInitConditions!(sdm::ShootingDataManager, icv::AbstractVector)
    resize!(sdm.ics, length(icv))
    sdm.ics .= icv 
    return nothing 
end

function SetFinConditions!(sdm::ShootingDataManager, fcv::AbstractVector)
    resize!(sdm.fcs, length(fcv))
    sdm.fcs .= fcv 
    return nothing
end

function SetPreallocatedVecs!(sdm::ShootingDataManager, size::Int)
    resize!(sdm.inVec, size)
    resize!(sdm.outVec, size)
    return nothing
end

function SetPreallocatedVecs!(sdm::ShootingDataManager, inVec::AbstractVector)
    n = length(inVec)
    resize!(sdm.inVec, n)
    resize!(sdm.outVec, n)
    sdm.inVec .= inVec 
    return nothing 
end

function GetInitGuessData(sdm::ShootingDataManager)
    return sdm.initGuessData
end

function GetInitConditions(sdm::ShootingDataManager)
    return sdm.ics
end

function GetFinConditions(sdm::ShootingDataManager)
    return sdm.fcs
end

function GetInputVec(sdm::ShootingDataManager)
    return sdm.inVec
end

function GetOutputVec(sdm::ShootingDataManager)
    return sdm.outVec 
end

