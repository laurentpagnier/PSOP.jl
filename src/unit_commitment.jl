export run_uc

function run_uc(
    ps::PSdata;
    shed_penalty = 10000,
    cur_penalty = 100,
    opt_tol = 0.01,
    node_limit = 500,
    load_share_res = 0.03,
    wind_share_res = 0.05,
    verbose = false,
    Nseg = 3,
)
    # create a gurobi model and set a few attibute
    model = Model(() -> Gurobi.Optimizer(GRB_ENV[]))
    if(!verbose)
        set_silent(model)
    end
    set_optimizer_attribute(model, "OptimalityTol", opt_tol)
    set_optimizer_attribute(model, "NodeLimit", node_limit)

    # define variables
    @variable(model, th[i=1:ps.Nbus,t=1:ps.Nt]) # phase angles
    @variable(model, ison[g=1:ps.Ngen,t=1:ps.Nt], Bin) 
    @variable(model, 0 <= startup[g=1:ps.Ngen,t=1:ps.Nt] <= 1)
    @variable(model, 0 <= shutdown[g=1:ps.Ngen,t=1:ps.Nt] <= 1)
    @variable(model, 0 <= gen[g=1:ps.Ngen,t=1:ps.Nt])
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

    cgmax = @constraint(model, gen + sync_res .<= ps.max_gen.*ison)
    cgmin = @constraint(model, ps.min_gen.*ison .<= gen)
    
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
    r = ps.ramping_rate
    pmax = ps.max_gen
    ot = ps.min_on_time
    dt = ps.min_down_time
    @constraint(model, ison[:,1] .== startup[:,1] - shutdown[:,1])
    for t=2:ps.Nt
        @constraint(model, ison[:,t] .== ison[:,t-1] + startup[:,t] - shutdown[:,t])
        # Ramping constraints
        @constraint(model, gen[:,t] - gen[:,t-1] .<= r.*ison[:,t-1] + pmax.*startup[:,t])
        @constraint(model, gen[:,t-1] - gen[:,t] .<= r.*ison[:,t-1] + pmax.*shutdown[:,t])
    end
    for t=1:ps.Nt
        for g=1:ps.Ngen
            @constraint(model, sum(startup[g,max(t-ot[g],1):t]) <= ison[g,t])
            @constraint(model, sum(shutdown[g,max(t-dt[g],1):t]) <= 1-ison[g,t])
        end
    end

    # Generation cost, here we approximate the quadratic cost function
    # by a piecewise linear function
    m, h = piecewise_cost(ps, Nseg = Nseg)
    for k=1:Nseg
        @constraint(model, gen_cost .>= m[:,k].*gen .+ h[:,k])
    end
    
    # Reserve constraints
    D = sum(ps.demand, dims=1)
    W = sum(ps.wind, dims=1)
    # here generator should be able the reserve in less than 10min, hence the 1/6.
    @constraint(model, sync_res .<= ps.ramping_rate / 6 .* ison)
    @constraint(model, req_res .<= sum(sync_res, dims=1)) # sum over generators
    @constraint(model, load_share_res*D + wind_share_res*W .<= sum(sync_res, dims=1))
    
    @objective(model, Min,
        sum(gen_cost) + shed_penalty * sum(shed) + cur_penalty * sum(cur)
        + sum(ps.startup_cost .* startup) + sum(ps.shutdown_cost .* shutdown)
    )
    
    optimize!(model)

    #fixing the decision variables and resolving
    # to obtain the LMPs
    temp1 = value.(ison)
    temp2 = value.(startup)
    temp3 = value.(shutdown)
    unset_binary.(ison)
    fix.(ison, temp1)
    fix.(startup, temp2, force=true)
    fix.(shutdown, temp3, force=true)
    optimize!(model)

    #lmp = dual.(cb)
    #mu = dual.(cfmax)
    #nu = dual.(cfmin)
    #lambda = [dual(cb[t]) for t=1:ps.Nt]

    #ison = value.(ison)
    #gen = value.(gen)
    #th = value.(th)
    #shed = value.(shed)
    #cur = value.(cur)
    return PSresult(value.(gen), value.(th), value.(shed), value.(cur),
        value.(ison), value.(startup), value.(shutdown), dual.(cb))
end


function piecewise_cost(ps::PSdata; Nseg::Int64 = 3)
    # approximate the quadratic cost by a piecewise linear function.
    # This function by designed usually overestimimates the cost.
    pmin = ps.min_gen
    pmax = ps.max_gen
    c0 = ps.on_cost
    c = ps.lin_cost
    q = ps.quad_cost
    dp = (pmax - pmin) / Nseg

    m = zeros(ps.Ngen, Nseg)
    h = zeros(ps.Ngen, Nseg)
    for g=1:ps.Ngen
        for k=1:Nseg
            m[g,k] = c[g] + q[g] * (2*pmin[g] + (2*k-1)*dp[g])
            h[g,k] = c0[g] - q[g] * (pmin[g] + (k-1)*dp[g])*(pmin[g] + k*dp[g])
        end
    end
    return m, h
end
