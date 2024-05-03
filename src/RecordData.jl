module RecordData
using CSV, DataFrames, ..Structs, GeoDataFrames, GeoFormatTypes, Arrow, ArchGDAL



function count_specific_land_use(agent::Agent, land_use_name::String)
    # return sum(parcel -> parcel.land_use_name == land_use_name, agent.land_parcels)
    return sum([p.area_ha for p in agent.land_parcels if p.land_use_name == land_use_name])
end


function record_initial_cap_dist(agents,land_uses,outputs_folder,suffix)
    println("recording cap dist")
    cap_dist = DataFrame(:uid=>[],:land_use_name=>[],:production_capability=>[],:active=>[])
    for agent in agents
        for land_use in land_uses
            pc = agent.plant_and_animal_capability_by_land_use[land_use.name] ^ 0.5 * agent.people_capability ^ 0.25 * agent.business_capability ^ 0.25
            if any([lp.land_use_name == land_use.name for lp in agent.land_parcels])
                push!(cap_dist, (agent.uid, land_use.name, pc, true))
            else
                push!(cap_dist, (agent.uid, land_use.name, pc, false))
            end
        end
    end
    Arrow.write(outputs_folder*"/cap_dist_"*string(suffix)*".arrow", cap_dist)
end

function export_agents(agents::Vector{Agent},land_uses::Vector{LandUse},land_info::DataFrame,timestep,outputs_folder)
    println("Timestep: ", timestep)
    # if (timestep==0) | (timestep==25)
    #     record_initial_cap_dist(agents,land_uses,outputs_folder,timestep)
    # end
    # for land_use in land_uses
    #     count_lu = sum(map(agent -> count_specific_land_use(agent, land_use.name), agents))
    #     # int_lu = sum(x->x.land_parcel.intensity,filter(x->x.land_parcel.land_use_name==land_use.name,agents))/count_lu
    #     print("Land use: ", land_use.name, " Count: ", count_lu, " \n")#, " Mitigation: ", mit_lu, " Profitability: ", prof_lu, " Environmental Damage Mitigated: ", env_lu, "\n")
    # end
    agent_vals = DataFrame(agents)
    agent_vals[:,:num_neighbours] = [length(a.neighbours) for a in agents]
    land_vals = DataFrame(lp for a in agents for lp in a.land_parcels)
    parcel_ids = land_vals[:,:uid]
    land_vals[:,:uid] = [a.uid for a in agents for lp in a.land_parcels]
    land_vals[:,:parcelid] = parcel_ids
    
    agent_vals = select!(agent_vals, Not(:land_parcels,:neighbours))
    
    df = outerjoin(agent_vals, land_vals, on=:uid)
    df = select(df, Not(:utility_barrier_matrix))
    for field in names(df)
        if occursin("by_land_use",field) # Dict{String,Float64} <: eltype(df[:,field])
            
            df[:,replace(field,"_by_land_use"=>"")] = [df[i,field][lu] for (i,lu) in enumerate(df[:,:land_use_name])]
            df = select!(df, Not(field))
        end
    end
    Arrow.write(outputs_folder*"/agents_"*string(timestep)*".arrow", df)
end

function export_land_use(agents::Vector{Agent},land_info::DataFrame,outputs_folder)
    land_vals = DataFrame(lp for a in agents for lp in a.land_parcels)
    with_geom = leftjoin(land_vals, land_info, on=:uid)
    with_geom[:,:x] = ArchGDAL.getx.(centroid.(with_geom[:,:geom]),Ref(0))
    with_geom[:,:y] = ArchGDAL.gety.(centroid.(with_geom[:,:geom]),Ref(0))
    only_lu_geom = select!(with_geom, [:land_use_name,:x,:y])
    Arrow.write(outputs_folder*"/land_use_init.arrow", only_lu_geom)
end

function export_extension_spend!(external_resources::Vector{ExternalResource},timestep,output_folder)
    external_resources_vals = DataFrame(external_resources)
    Arrow.write(output_folder*"/extension_"*string(timestep)*".arrow", external_resources_vals)
    for er in external_resources
        er.run_this_timestep = false
    end
end


function export_catchments(catchments::Array{Catchment,1},timestep,output_folder)
    catchment_vals = DataFrame(catchments)
    select!(catchment_vals, Not(:land_parcels))
    Arrow.write(output_folder*"/catchments_"*string(timestep)*".arrow", catchment_vals)
end

end