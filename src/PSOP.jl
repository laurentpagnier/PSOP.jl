module PSOP

using Gurobi
using JuMP
using SparseArrays
using LinearAlgebra
using HDF5

struct PSdata
    gen_loc::Vector{Int64} # where (i.e. bus ids) the generators are located
    wind_loc::Vector{Int64}
    min_gen::Vector{Float64} # all values in this structure are stored in pu
    max_gen::Vector{Float64}
    line_id::Matrix{Int64}
    line_susceptance::Vector{Float64}
    line_limit::Vector{Float64}
    demand::Matrix{Float64}
    wind::Matrix{Float64}
    ramping_rate::Vector{Float64} # in pu/hour
    lin_cost::Vector{Float64}
    quad_cost::Vector{Float64}
    on_cost::Vector{Float64} # cost for being on duty
    min_on_time::Vector{Int64}
    min_down_time::Vector{Int64}
    Nbus::Int64
    Nline::Int64
    Ngen::Int64
    Nwind::Int64
    Nt::Int64
    sb::Float64 # base MVA
end

include("data_handler.jl")
include("optimal_power_flow.jl")
include("unit_commitment.jl")
include("real_time_dispatch.jl")

end
