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
)
    model = Model(Gurobi.Optimizer)
    if(!verbose)
        set_silent(model)
    end
    set_optimizer_attribute(model, "OptimalityTol", opt_tol)
    set_optimizer_attribute(model, "NodeLimit", node_limit)

    # define variables
    @variable(model, th[i=1:ps.Nbus,t=1:ps.Nt])
    @variable(model, ison[g=1:ps.Ngen,t=1:ps.Nt], Bin)
    @variable(model, 0 <= startup[g=1:ps.Ngen,t=1:ps.Nt] <= 1)
    @variable(model, 0 <= shutdown[g=1:ps.Ngen,t=1:ps.Nt] <= 1)
    @variable(model, 0 <= gen[g=1:ps.Ngen,t=1:ps.Nt])
    @variable(model, 0 <= gen_cost[g=1:ps.Ngen,t=1:ps.Nt])
    @variable(model, 0 <= spin_res[g=1:ps.Ngen,t=1:ps.Nt])
    @variable(model, 0 <= shed[i=1:ps.Nbus,t=1:ps.Nt] <= ps.demand[i,t]) # load shedding
    @variable(model, 0 <= cur[w=1:ps.Nwind,t=1:ps.Nt] <= ps.wind[w,t]) # wind curtailment

    # slack ref constraint
    for t=1:ps.Nt
        @constraint(model, th[1,t] == 0)
    end

    for t=1:ps.Nt
        for g=1:ps.Ngen
            c1 = @constraint(model, gen[g,t] + spin_res[g,t] <= ps.max_gen[g]*ison[g,t])
            c2 = @constraint(model, ps.min_gen[g]*ison[g,t] <= gen[g,t])
        end
    end

    # line limit constraints
    line_cons = []
    b = ps.line_susceptance
    id1 = ps.line_id[:,1]
    id2 = ps.line_id[:,2]
    for t=1:ps.Nt
        for k=1:ps.Nline
            if ps.line_limit[k] > 0
                c1 = @constraint(model, b[k] * (th[id1[k],t] - th[id2[k],t]) <= ps.line_limit[k])
                c2 = @constraint(model, b[k] * (th[id2[k],t] - th[id1[k],t]) <= ps.line_limit[k])
                push!(line_cons, c1)
                push!(line_cons, c2)
            end
        end
    end
    
    # power injection constraints
    inc = sparse([ps.line_id[:,1]; ps.line_id[:,2]], [1:ps.Nline; 1:ps.Nline],
        [-ones(ps.Nline); ones(ps.Nline)])
    B = inc * (ps.line_susceptance .* inc')
    gen2bus = sparse(ps.gen_loc, 1:ps.Ngen, ones(ps.Ngen), ps.Nbus, ps.Ngen)
    wind2bus = sparse(ps.wind_loc, 1:ps.Nwind, ones(ps.Nwind), ps.Nbus, ps.Nwind)
    balance_cons = []
    for t=1:ps.Nt
        c = @constraint(model, gen2bus * gen[:,t] - B * th[:,t] .== ps.demand[:,t]
            - shed[:,t] - wind2bus * (ps.wind[:,t] - cur[:,t]))
        push!(balance_cons, c)
    end
    
    
    ot = ps.min_on_time
    dt = ps.min_down_time
    for t=1:ps.Nt
        for g=1:ps.Ngen
            @constraint(model, sum(startup[g,max(t-ot[g],1):t])  <= ison[g,t])
            @constraint(model, sum(shutdown[g,max(t-dt[g],1):t])  <= 1-ison[g,t])
        end
    end
    
    for g=1:ps.Ngen
        @constraint(model, ison[g,1] == startup[g,1] - shutdown[g,1])
    end
    for t=2:ps.Nt
        for g=1:ps.Ngen
            r = ps.ramping_rate[g]
            pmax = ps.max_gen[g]
            @constraint(model, ison[g,t] == ison[g,t-1] + startup[g,t] - shutdown[g,t])
            # Ramping constraints
            @constraint(model, gen[g,t] - gen[g,t-1] <= r*ison[g,t-1] + pmax*startup[g,t])
            @constraint(model, gen[g,t-1] - gen[g,t] <= r*ison[g,t-1] + pmax*shutdown[g,t])
        end
    end

    # Generation cost
    m, h = piecewise_cost(ps)
    for t=1:ps.Nt
        for g=1:ps.Ngen
            for k=1:size(m,2)
                @constraint(model, gen_cost[g,t] >= m[g,k]*gen[g,t] + h[g,k])
            end
        end
    end
    
    # Reserve constraints
    D = sum(ps.demand, dims=1)
    W = sum(ps.wind, dims=1)
    for t=1:ps.Nt
        for g=1:ps.Ngen
            r = ps.ramping_rate[g]
            @constraint(model, spin_res[g,t] <= r/6 * ison[g,t])
        end
    end
    for t=1:ps.Nt
        @constraint(model, load_share_res*D[t] + wind_share_res*W[t] <= sum(spin_res[:,t]))
    end
    
    @objective(model, Min,
        sum(gen_cost)
        + shed_penalty * sum(shed)
        + cur_penalty * sum(cur)
        + sum(ps.on_cost .* ison))
    
    optimize!(model)

    #fixing the decision variables and resolving
    # to obtain the LMPs
    temp1 = value.(ison)
    temp2 = value.(startup)
    temp3 = value.(shutdown)
    unset_binary.(ison)
    #unset_binary.(startup)
    #unset_binary.(shutdown)
    fix.(ison, temp1)
    fix.(startup, temp2, force=true)
    fix.(shutdown, temp3, force=true)
    optimize!(model)

    lmp = [dual.(balance_cons[t][i]) for i=1:ps.Nbus, t=1:1:ps.Nt]

    ison = value.(ison)
    gen = value.(gen)
    th = value.(th)
    shed = value.(shed)
    cur = value.(cur)
    return gen, th, shed, cur, ison, lmp
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
            h[g,k] = c0[g] - q[g] * (pmin[g] + (k-1)*dp[g])*(pmin[g]+k*dp[g])
        end
    end
    return m, h
end
