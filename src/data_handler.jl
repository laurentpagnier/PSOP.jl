export import_data, save_as_hdf5, import_data_from_hdf5

function import_data(
    dirname::String
)
    # load the structure of the power system
    text = read(dirname*"case.m", String)
    mva = collect.(findall("mpc.baseMVA", text))[1][end]+1
    b = collect.(findall(r"mpc.bus( |=)", text))[1][end]+1
    l = collect.(findall("mpc.branch", text))[1][end]+1
    g = collect.(findall(r"mpc.gen( |=)", text))[1][end]+1
    gc = collect.(findall("mpc.gencost", text))[1][end]+1
    gc = collect.(findall("mpc.gencost", text))[1][end]+1
    bs = collect.(findall("]", text)) .|> (x) -> x[1]
    rs = collect.(findall("\n", text)) .|> (x) -> x[1]-1
    
    mva_data = text[mva:rs[findfirst(rs .> mva)]]    
    bus_data = text[b:bs[findfirst(bs .> b)]]
    line_data = text[l:bs[findfirst(bs .> l)]]
    gen_data = text[g:bs[findfirst(bs .> g)]]
    gencost_data = text[gc:bs[findfirst(bs .> gc)]]
    
    eval(Meta.parse("basemva $(mva_data)")) # the equal sign is included in the second term
    eval(Meta.parse("bus $bus_data"))
    eval(Meta.parse("line $line_data"))
    eval(Meta.parse("gen $gen_data"))
    eval(Meta.parse("gencost $gencost_data"))

    # load the demand profile
    text = read(dirname*"load_profile.csv", String)
    text = replace(text, ","=>"\t")
    eval(Meta.parse("demand = [$text]"))
    
    # load the wind profile
    text = read(dirname*"wind_profile.csv", String)
    text = replace(text, ","=>"\t")
    eval(Meta.parse("wind = [$text]"))

    wind_loc = collect(1:size(wind,1)) # !!!!!!!!!! quick fixe, NEED TO BE CHANGED 
    min_on_time = 5*ones(size(gen,1)) # !!!!!!!!!! quick fixe, NEED TO BE CHANGED 
    min_down_time = 5*ones(size(gen,1)) # !!!!!!!!!! quick fixe, NEED TO BE CHANGED 
    
    return PSdata(Int64.(gen[:,1]), wind_loc, gen[:,10] / basemva, gen[:,9] / basemva,
        Int64.(line[:,1:2]), 1 ./ line[:,4], line[:,6] / basemva, demand / basemva,
        wind / basemva, gencost[:,6], gencost[:,5], min_on_time, min_down_time,
        size(bus,1), size(line,1), size(gen,1), size(wind,1), size(demand,2), basemva)
end


function import_data_from_hdf5(
    filename::String
)
    # This function is not foolproof, it expects all the fields to
    # be provided.
    raw_data = h5read(filename, "/")
    fields = fieldnames(PSdata)
    data = Tuple(raw_data[String(f)] for f in fields) 
    return PSdata(data...)
end


function save_as_hdf5(
    filename::String,
    ps::PSdata,
)
    fields = fieldnames(PSdata)
    fid = h5open(filename, "w")
    for f in fields
        fid[String(f)] = getproperty(ps, f)
    end
    close(fid)
end


