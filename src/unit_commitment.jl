export run_uc

function run_uc(ps::PSdata;
    shed_penalty = 10000,
    cur_penalty = 1000,
)
    model = Model(Gurobi.Optimizer)
    set_silent(model)

    # define variables
    @variable(model, th[i=1:ps.Nbus,t=1:ps.Nt])
    @variable(model, ison[g=1:ps.Ngen,t=1:ps.Nt], Bin)
    @variable(model, startup[g=1:ps.Ngen,t=1:ps.Nt], Bin)
    @variable(model, shutdown[g=1:ps.Ngen,t=1:ps.Nt], Bin)
    @variable(model, gen[g=1:ps.Ngen,t=1:ps.Nt])
    @variable(model, 0 <= shed[i=1:ps.Nbus,t=1:ps.Nt] <= ps.demand[i,t]) # load shedding
    @variable(model, 0 <= cur[w=1:ps.Nwind,t=1:ps.Nt] <= ps.demand[w,t]) # wind curtailment

    # slack ref constraint
    for t=1:ps.Nt
        @constraint(model, th[1,t] == 0)
    end

    for t=1:ps.Nt
        for g=1:ps.Ngen
            c1 = @constraint(model, gen[g,t] <= ps.max_gen[g]*ison[g,t])
            c2 = @constraint(model, ps.min_gen[g]*ison[g,t] <= gen[g,t])
        end
    end

    # line limit constraints
    line_cons = []
    for t=1:ps.Nt
        for k=1:ps.Nline
            if ps.line_limit[k] > 0
                c1 = @constraint(model, ps.line_susceptance[k] *
                    (th[ps.line_id[k,1],t] - th[ps.line_id[k,2],t]) <=  ps.line_limit[k])
                c2 = @constraint(model, ps.line_susceptance[k] *
                    (th[ps.line_id[k,2],t] - th[ps.line_id[k,1],t]) <=  ps.line_limit[k])
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
        c = @constraint(model, gen2bus * gen[:,t] - B * th[:,t] .== ps.demand[:,t]
            - shed[:,t] - wind2bus*(ps.wind[:,t]-cur[:,t]))
        push!(balance_cons, c)
    end
    
    
    for t=1:ps.Nt
        for g=1:ps.Ngen
            @constraint(model, sum(startup[g,max(t-ps.min_on_time[g],1):t])  <= ison[g,t])
            @constraint(model, sum(shutdown[g,max(t-ps.min_down_time[g],1):t])  <= 1-ison[g,t])
        end
    end
    

    for t=2:ps.Nt
        for g=1:ps.Ngen
            @constraint(model, ison[g,t] == ison[g,t-1] + startup[g,t] - shutdown[g,t])
        end
    end
    
    @objective(model, Min, sum(ps.lin_cost.* gen)
        + sum(ps.quad_cost .* gen.^2)
        + shed_penalty*sum(shed))
    
    optimize!(model)

    ison = value.(ison)
    gen = value.(gen)
    th = value.(th)
    shed = value.(shed)
    cur = value.(cur)
    return gen, th, shed, cur
end
