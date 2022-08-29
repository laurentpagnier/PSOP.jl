export run_std_dc_opf, run_ptdf_dc_opf


function run_std_dc_opf(
    ps::PSdata;
    t = 1
)
    model = Model(() -> Gurobi.Optimizer(GRB_ENV))
    set_silent(model)
    
    # define variables
    @variable(model, ps.min_gen[g] <= gen[g=1:ps.Ngen] <= ps.max_gen[g] )
    @variable(model, th[i=1:ps.Nbus])

    # slack ref constraint
    @constraint(model, th[1] == 0)

    # line limit constraints
    line_cons = [];
    for k=1:ps.Nline
        if ps.line_limit[k] > 0
            c1 = @constraint(model, ps.line_susceptance[k] *
                (th[ps.line_id[k,1]] - th[ps.line_id[k,2]]) <=  ps.line_limit[k])
            c2 = @constraint(model, ps.line_susceptance[k] *
                (th[ps.line_id[k,2]] - th[ps.line_id[k,1]]) <=  ps.line_limit[k])
            push!(line_cons,c1)
            push!(line_cons,c2)
        end
    end
    
    # power injection constraints
    inc = sparse([ps.line_id[:,1]; ps.line_id[:,2]], [1:ps.Nline; 1:ps.Nline],
        [-ones(ps.Nline); ones(ps.Nline)])
    B = inc * (ps.line_susceptance .* inc')
    gen2bus = sparse(ps.gen_loc, 1:ps.Ngen, ones(ps.Ngen), ps.Nbus, ps.Ngen)
    balance_cons = @constraint(model, gen2bus * gen - B * th .==  ps.demand[:,t])

    @objective(model, Min, sum(ps.lin_cost.* gen)
        + sum(ps.quad_cost .* gen.^2))

    optimize!(model)
    
    if(termination_status(model) == MOI.INFEASIBLE)
        println("The problem is infeasible")
        return nothing, nothing, nothing
    else
        lmp = dual.(balance_cons)
        return value.(th), value.(gen), lmp
    end
end


function run_ptdf_dc_opf(
    ps::PSdata;
    t = 1
)
    # This version of is there mostly for testing purposes, it is significantly
    # less efficient the std_dc_opf for large systems
    model = Model(() -> Gurobi.Optimizer(GRB_ENV))
    set_silent(model)
    
    ptdf = create_ptdf_matrix(ps)
    
    # define variables
    @variable(model, ps.min_gen[g] <= gen[g=1:ps.Ngen] <= ps.max_gen[g] )

    # balance constraint
    c = @constraint(model, sum(gen) == sum(ps.demand[:,t]))
    
    # line limit constraints
    gen2bus = sparse(ps.gen_loc, 1:ps.Ngen, ones(ps.Ngen), ps.Nbus, ps.Ngen)
    id_lim = findall(ps.line_limit .> 0)
    if !isempty(id_lim)
        c_line1 = @constraint(model, ptdf[id_lim,:] * (gen2bus * gen 
            - ps.demand[:,t]) .<= ps.line_limit[id_lim])
        c_line2 = @constraint(model, ptdf[id_lim,:] * (gen2bus * gen 
            - ps.demand[:,t]) .>= -ps.line_limit[id_lim])
    else
        c_line1 = []
        c_line2 = []
    end
    @objective(model, Min, sum(ps.lin_cost.* gen)
        + sum(ps.quad_cost .* gen.^2))
    
    optimize!(model)
    if(termination_status(model) == MOI.INFEASIBLE)
        println("The problem is infeasible")
        return nothing, nothing
    else
        # here the price consist of 2 terms: the means price and
        # local deviations due to congestion
        lmp = dual(c) .+ ptdf[id_lim,:]' * (dual.(c_line1) + dual.(c_line2))
        return value.(gen), lmp
    end
end


function create_ptdf_matrix(ps::PSdata; tol = 1E-11)
    inc = sparse([ps.line_id[:,1]; ps.line_id[:,2]], [1:ps.Nline; 1:ps.Nline],
        [-ones(ps.Nline); ones(ps.Nline)])
    B = inc * (ps.line_susceptance .* inc')
    eig = eigen(Matrix{Float64}(B))
    d =  eig.values
    U = eig.vectors
    d[d .> tol] .= 1 ./ d[d .> tol]
    Bt = U * (d .* U')
    return ps.line_susceptance .* (inc' * Bt)
end
