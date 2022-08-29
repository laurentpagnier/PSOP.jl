export run_rt

function run_rt(
    ps::PSdata,
    scuc_results;
    shed_penalty = 1000,
    cur_penalty = 100,
    redispatch_penalty = 50,
)
    model = Model(() -> Gurobi.Optimizer(GRB_ENV))
    #set_silent(model)


    # extract scuc results    
    uc_gen, _ , _ , _ , uc_ison, _ = scuc_results
    rt_per_uc = Int64.(size(ps.demand,2) / size(uc_gen,2))

    # define variables
    @variable(model, th[i=1:ps.Nbus,t=1:ps.Nt])
    @variable(model, gen[g=1:ps.Ngen,t=1:ps.Nt])
    @variable(model, 0 <= shed[i=1:ps.Nbus,t=1:ps.Nt] <= ps.demand[i,t]) # load shedding
    @variable(model, 0 <= cur[w=1:ps.Nwind,t=1:ps.Nt] <= ps.wind[w,t]) # wind curtailment

    # slack ref constraint
    for t=1:ps.Nt
        @constraint(model, th[1,t] == 0)
    end

    for t=1:ps.Nt
        for g=1:ps.Ngen
            uc_t = Int64.(floor((t-1)/rt_per_uc) + 1)
            c1 = @constraint(model, gen[g,t] <= ps.max_gen[g]*uc_ison[g,uc_t] - uc_gen[g,uc_t])
            c2 = @constraint(model, ps.min_gen[g]*uc_ison[g,uc_t] - uc_gen[g,uc_t] <= gen[g,t])
        end
    end

    # line limit constraints
    line_cons = []
    for t=1:ps.Nt
        for k=1:ps.Nline
            if ps.line_limit[k] > 0
                c1 = @constraint(model, ps.line_susceptance[k] *
                    (th[ps.line_id[k,1],t] - th[ps.line_id[k,2],t]) <= ps.line_limit[k])
                c2 = @constraint(model, ps.line_susceptance[k] *
                    (th[ps.line_id[k,2],t] - th[ps.line_id[k,1],t]) <= ps.line_limit[k])
                push!(line_cons,c1)
                push!(line_cons,c2)
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
        uc_t = Int64.(floor((t-1)/rt_per_uc) + 1)
        c = @constraint(model, gen2bus * (gen[:,t] + uc_gen[:,uc_t]) 
            - B * th[:,t] .== ps.demand[:,t]
            - shed[:,t] - wind2bus*(ps.wind[:,t] - cur[:,t]))
        push!(balance_cons, c)
    end
    
    @objective(model, Min,
        sum(ps.lin_cost.* gen)
        + sum(ps.quad_cost .* gen.^2)
        + shed_penalty * sum(shed)
        + cur_penalty * sum(cur))
    
    optimize!(model)

    lmp = [dual.(balance_cons[t][i]) for i=1:ps.Nbus, t=1:1:ps.Nt]

    gen = value.(gen)
    th = value.(th)
    shed = value.(shed)
    cur = value.(cur)
    return gen, th, shed, cur, lmp
end
