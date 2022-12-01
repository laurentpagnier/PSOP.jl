export run_rt, create_real_time_from_uc

function run_rt(
    ps::PSdata,
    scuc_results::PSresult;
    shed_penalty = 1000,
    cur_penalty = 100,
    redispatch_penalty = 50,
    opt_tol = 0.01,
    node_limit = 500,
    verbose = false,
    quad_cost = false,
    Nseg = 3,
)
    # extract scuc results    
    uc_gen, _ , _ , _ , uc_ison, _ = scuc_results
    rt_per_uc = Int64.(ps.Nt / size(uc_gen,2))
    
    # extend scuc results to same time resolution as real time
    uc_gen = kron(scuc_results.gen, ones(1, rt_per_uc))
    uc_ison = kron(scuc_results.ison, ones(1, rt_per_uc))
    uc_shartup = kron(scuc_results.startup, ones(1, rt_per_uc))
    uc_shutdown = kron(scuc_results.shutdown, ones(1, rt_per_uc))
    
    # create a gurobi model and set a few attibute
    model = Model(() -> Gurobi.Optimizer(GRB_ENV[]))
    if(!verbose)
        set_silent(model)
    end
    set_optimizer_attribute(model, "OptimalityTol", opt_tol)
    set_optimizer_attribute(model, "NodeLimit", node_limit)


    # define variables
    @variable(model, th[i=1:ps.Nbus,t=1:ps.Nt]) # phase angles
    @variable(model, 0 <= pos_redispatch[g=1:ps.Ngen,t=1:ps.Nt])
    @variable(model, 0 <= neg_redispatch[g=1:ps.Ngen,t=1:ps.Nt])
    @variable(model, 0 <= gen_cost[g=1:ps.Ngen,t=1:ps.Nt]) # generation cost
    # reserve provided by synchronous (i.e. on duty) generators 
    @variable(model, 0 <= sync_res[g=1:ps.Ngen,t=1:ps.Nt]) 
    # reserve provided by wind farm
    #@variable(model, 0 <= wind_res[g=1:ps.Nwind,t=1:ps.Nt]) # 
    # reserve provided by fast synchronizing generators
    #@variable(model, 0 <= fast_res[g=1:ps.Ngen,t=1:ps.Nt]) # fast re
    @variable(model, 0 <= shed[i=1:ps.Nbus,t=1:ps.Nt] <= ps.demand[i,t]) # load shedding
    @variable(model, 0 <= cur[w=1:ps.Nwind,t=1:ps.Nt] <= ps.wind[w,t]) # wind curtailment

    # Counterintuitively, it seems that introducing additional intermediate
    # variables helps gurobi to solve the problem faster
    @variable(model, 0 <= gen[g=1:ps.Ngen,t=1:ps.Nt])
    @variable(model, flow[k=1:ps.Nline,t=1:ps.Nt])
    @variable(model, wind_gen[k=1:ps.Nwind,t=1:ps.Nt])
    @variable(model, demand[i=1:ps.Nbus,t=1:ps.Nt])
    @variable(model, 0 <= req_res[t=1:ps.Nt]) # required reserve

    # slack ref constraint
    # 1st bus is treated as the slack bus
    for t=1:ps.Nt
        @constraint(model, th[1,t] == 0)
    end

    cd = @constraint(model, demand .== ps.demand - shed)
    cw = @constraint(model, wind_gen .== ps.wind - cur)
    cg = @constraint(model, gen .== uc_gen + pos_redispatch - neg_redispatch)

    cgmax = @constraint(model, gen + sync_res .<= ps.max_gen.*uc_ison)
    cgmin = @constraint(model, ps.min_gen.*uc_ison .<= gen)
    
    line_cons = []
    b = ps.line_susceptance
    id1 = ps.line_id[:,1]
    id2 = ps.line_id[:,2]
    cf = @constraint(model, flow .== b .* (th[id1,:] - th[id2,:]))
    cfmax = @constraint(model, flow .<= ps.line_limit)
    cfmin = @constraint(model, flow .>= -ps.line_limit)
    
    # power balance constraints
    gen2bus = sparse(ps.gen_loc, 1:ps.Ngen, ones(ps.Ngen), ps.Nbus, ps.Ngen)
    wind2bus = sparse(ps.wind_loc, 1:ps.Nwind, ones(ps.Nwind), ps.Nbus, ps.Nwind)
    from_bus = sparse(id1, 1:ps.Nline, ones(ps.Nline), ps.Nbus, ps.Nline)
    to_bus = sparse(id2, 1:ps.Nline, ones(ps.Nline), ps.Nbus, ps.Nline)
    # local balance
    cb = @constraint(model, gen2bus * gen + wind2bus * wind_gen
            + to_bus * flow - from_bus * flow .== demand)
    # global balance (optional) might help gurobi
    #cb2 = @constraint(model, sum(gen, dims=1)
    #    + sum(wind_gen, dims=1) .== sum(demand, dims=1))

    # decision variable and ramping constraints
    r = ps.ramping_rate / rt_per_uc # here UC is assumed to be hourly
    pmax = ps.max_gen
    ot = ps.min_on_time
    dt = ps.min_down_time
    @constraint(model, ison[:,1] .== startup[:,1] - shutdown[:,1])
    for t=2:ps.Nt
        #@constraint(model, ison[:,t] .== ison[:,t-1] + startup[:,t] - shutdown[:,t])
        # Ramping constraints
        @constraint(model, gen[:,t] - gen[:,t-1] .<= r.*uc_ison[:,t-1] + pmax.*startup[:,t])
        @constraint(model, gen[:,t-1] - gen[:,t] .<= r.*uc_ison[:,t-1] + pmax.*shutdown[:,t])
    end
    #=
    for t=1:ps.Nt
        for g=1:ps.Ngen
            @constraint(model, sum(startup[g,max(t-ot[g],1):t]) <= ison[g,t])
            @constraint(model, sum(shutdown[g,max(t-dt[g],1):t]) <= 1-ison[g,t])
        end
    end
    =#
    # Generation cost, here we approximate the quadratic cost function
    # by a piecewise linear function
    m, h = PSOP.piecewise_cost(ps, Nseg = Nseg)
    for k=1:Nseg
        @constraint(model, gen_cost .>= m[:,k].*gen .+ h[:,k])
    end
    
    # Reserve constraints
    #D = sum(ps.demand, dims=1)
    #W = sum(ps.wind, dims=1)
    # here generator should be able the reserve in less than 10min, hence the 1/6.
    #@constraint(model, sync_res .<= ps.ramping_rate / 6 .* ison)
    #@constraint(model, req_res .<= sum(sync_res, dims=1)) # sum over generators
    #@constraint(model, load_share_res*D + wind_share_res*W .<= sum(sync_res, dims=1))
    
    @objective(model, Min,
        sum(gen_cost) + redispatch_penalty * (sum(pos_redispatch) + sum(neg_redispatch)) 
        + shed_penalty * sum(shed) + cur_penalty * sum(cur)
    )
    
    optimize!(model)


    lmp = dual.(cb)
    mu = dual.(cfmax)
    nu = dual.(cfmin)
    #lambda = [dual(cb[t]) for t=1:ps.Nt]

    #ison = value.(ison)
    gen = value.(gen)
    th = value.(th)
    shed = value.(shed)
    return gen, th, shed, cur, ison, lmp
end


function create_real_time_from_uc(
    ps::PSdata,
    demand::Matrix{Float64},
    wind::Matrix{Float64},
)    
    rt_ps = PSOP.PSdata(ps.gen_loc, ps.wind_loc, ps.min_gen, ps.max_gen, ps.line_id, ps.line_susceptance,
    ps.line_limit, demand, wind, ps.ramping_rate, ps.lin_cost, ps.quad_cost, ps.on_cost, ps.startup_cost,
    ps.shutdown_cost, ps.min_on_time, ps.min_down_time, ps.Nbus, ps.Nline, ps.Ngen, ps.Nwind, size(wind,2),
    ps.sb)
end
