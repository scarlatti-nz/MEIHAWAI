module Population

using ..Structs, Distributions, DataFrames, GeoDataFrames, XLSX, ArchGDAL



function get_land_attributes(area_ha,allowed_land_uses,land_use_inputs,dairy_typology_inputs,row)
    #   Initalise dictionaries for parameter storage
    max_revenue_by_land_use = Dict{String,Float64}()
    min_intensity_by_land_use = Dict{String,Float64}()
    max_intensity_by_land_use = Dict{String,Float64}()
    n_mitigation_potential_by_land_use = Dict{String,Float64}()
    p_mitigation_potential_by_land_use = Dict{String,Float64}()
    nitrogen_loss_coefficients_by_land_use = Dict{String,Tuple{Float64,Float64,Float64}}()
    phosphorous_loss_coefficients_by_land_use = Dict{String,Tuple{Float64,Float64,Float64}}()
    nitrous_oxide_loss_coefficients_by_land_use = Dict{String,Tuple{Float64,Float64,Float64}}()
    sediment_intensity_coefficients_by_land_use = Dict{String,Tuple{Float64,Float64}}()
    sediment_loss_potenial_by_land_use = Dict{String,Float64}()
    carbon_sequestration_adjustment_by_land_use = Dict{String,Float64}()
    typecode = row["dairy_typology_typecode"]
    if !(typecode in dairy_typology_inputs[:,:typecode])
        println("Warning: ",row[:uid]," has typecode ",typecode," which is not in dairy_typology")
        typecode = 1
    end
    #   Fill the dictionaries with parameter values for each land use, different land uses as the keys
    for land_use in allowed_land_uses
        #   Map values for revenue to the land use keys, divide by area to get average from sum
        max_revenue_by_land_use[land_use] = row["economic_indicator_revenue-"*land_use*"_value"] / area_ha
        min_intensity_by_land_use[land_use] = dairy_typology_inputs[(dairy_typology_inputs[:,:typecode] .== typecode) .&& (dairy_typology_inputs[:,"Land use"] .== land_use),:Intensity_lower_limit][1]
        max_intensity_by_land_use[land_use] = dairy_typology_inputs[(dairy_typology_inputs[:,:typecode] .== typecode) .&& (dairy_typology_inputs[:,"Land use"] .== land_use),:Intensity_upper_limit][1]
        n_mitigation_potential_by_land_use[land_use] = dairy_typology_inputs[(dairy_typology_inputs[:,:typecode] .== typecode) .&& (dairy_typology_inputs[:,"Land use"] .== land_use),:N_mitigation_potential][1]
        p_mitigation_potential_by_land_use[land_use] = dairy_typology_inputs[(dairy_typology_inputs[:,:typecode] .== typecode) .&& (dairy_typology_inputs[:,"Land use"] .== land_use),:P_mitigation_potential][1]
        #   Map values for other keys
        nitrogen_loss_coefficients_by_land_use[land_use] = (dairy_typology_inputs[(dairy_typology_inputs[:,:typecode] .== typecode) .&& (dairy_typology_inputs[:,"Land use"] .== land_use),:Constant_N][1],
                                                            dairy_typology_inputs[(dairy_typology_inputs[:,:typecode] .== typecode) .&& (dairy_typology_inputs[:,"Land use"] .== land_use),:Linear_N][1],
                                                            dairy_typology_inputs[(dairy_typology_inputs[:,:typecode] .== typecode) .&& (dairy_typology_inputs[:,"Land use"] .== land_use),:Quadratic_N][1])
        phosphorous_loss_coefficients_by_land_use[land_use] = (dairy_typology_inputs[(dairy_typology_inputs[:,:typecode] .== typecode) .&& (dairy_typology_inputs[:,"Land use"] .== land_use),:Constant_P][1],
                                                                dairy_typology_inputs[(dairy_typology_inputs[:,:typecode] .== typecode) .&& (dairy_typology_inputs[:,"Land use"] .== land_use),:Linear_P][1],
                                                                dairy_typology_inputs[(dairy_typology_inputs[:,:typecode] .== typecode) .&& (dairy_typology_inputs[:,"Land use"] .== land_use),:Quadratic_P][1])
        nitrous_oxide_loss_coefficients_by_land_use[land_use] = (dairy_typology_inputs[(dairy_typology_inputs[:,:typecode] .== typecode) .&& (dairy_typology_inputs[:,"Land use"] .== land_use),:Constant_N2O][1],
                                                                dairy_typology_inputs[(dairy_typology_inputs[:,:typecode] .== typecode) .&& (dairy_typology_inputs[:,"Land use"] .== land_use),:Linear_N2O][1],
                                                                dairy_typology_inputs[(dairy_typology_inputs[:,:typecode] .== typecode) .&& (dairy_typology_inputs[:,"Land use"] .== land_use),:Quadratic_N2O][1])
        sediment_intensity_coefficients_by_land_use[land_use] = (dairy_typology_inputs[(dairy_typology_inputs[:,:typecode] .== typecode) .&& (dairy_typology_inputs[:,"Land use"] .== land_use),:Constant_Sediment][1],
                                                                dairy_typology_inputs[(dairy_typology_inputs[:,:typecode] .== typecode) .&& (dairy_typology_inputs[:,"Land use"] .== land_use),:Linear_Sediment][1])
        # construct the sediment column name from the coverclass in the inputs spreadsheet
        sediment_loss_potenial_by_land_use[land_use] = row["sediment_load_" * land_use_inputs[land_use_inputs[:,:Crop] .== land_use,:Sediment_cover_class][1]] / area_ha
        # scale carbon sequestration by 45000 to convert from revenue back to yield
        carbon_sequestration_adjustment_by_land_use[land_use] = row["economic_indicator_revenue-production-forestry_value"] / (area_ha * 45000)
    end

    #   Return dictionaries with parameters mapped to different land uses
    return typecode, min_intensity_by_land_use, max_intensity_by_land_use, max_revenue_by_land_use, n_mitigation_potential_by_land_use, p_mitigation_potential_by_land_use, nitrogen_loss_coefficients_by_land_use, phosphorous_loss_coefficients_by_land_use, nitrous_oxide_loss_coefficients_by_land_use, sediment_intensity_coefficients_by_land_use, sediment_loss_potenial_by_land_use, carbon_sequestration_adjustment_by_land_use
end

function read_cost_to_switch_matrix(name,cost_to_switch)

    #   Read in the cost switching matrix - cost to switch land use for all possible combindations
    

    #   Filter rows where the first column value matches input
    values = filter(row -> row[:to] == name, cost_to_switch)
    
    #   map column name from cost_to_switch to first row in values, and store dictionary of values
    cost_to_switch_dict = Dict(k=>v for (k,v) in zip(names(cost_to_switch)[2:end],values[1,2:end]))

    return cost_to_switch_dict
end

function read_spillover_learning_matrix(name,spillover_learning_matrix)
    #   Read in transferability dataframe
    
    
    #   Filter rows where the first column value matches input
    values = filter(row -> row[:from] == name, spillover_learning_matrix)
    
    #   map column name from spill_over to first row in values, and store dictionary of values
    spillover_learning_dict = Dict(k=>v for (k,v) in zip(names(spillover_learning_matrix)[2:end],values[1,2:end]))
    return spillover_learning_dict
end

function read_carbon_tables(name)
    #   If 'Forestry' is not found in name, return an array of 100 zeros, if not forestry related, no carbon_sequestration
    if ! occursin("forestry",name)
        return zeros(100)
    end
    carbon_table = DataFrame(XLSX.readtable("data/Model inputs.xlsx","carbon_sequestration"))
    values = carbon_table[:,name]
    return values
end


function initialise_land_uses(;enable_additional_uses=[])

    #   Read in datafile and store in data frame
    land_use_inputs = DataFrame(XLSX.readtable("data/Model inputs.xlsx","land_use_parameters"))
 
    #   Check if additional_uses is empty, if yes, only keep rows where active column is yes
    if isempty(enable_additional_uses)
        filter!(row -> row[:Active] == "Yes", land_use_inputs)
    #   If additional_uses, retain rows where 'Crop' column is in additional uses
    else
        filter!(row -> row[:Crop] in enable_additional_uses, land_use_inputs)
    end
    
    #   Land_use_names is the value from 'Crop' column
    land_use_names = land_use_inputs[:,:Crop]
    land_use_categories = land_use_inputs[:,:Category]

    #   Retrieves relevant numerical data from specified column and converts to float
    fixed_costs = float.(land_use_inputs[:,"Fixed Expenses (\$/ha)"])
    years_to_yield = float.(land_use_inputs[:,"Years to yield (avg)"])
    crop_lifetime = float.(land_use_inputs[:,"Lifetime"])
    revenue_adjustment = float.(land_use_inputs[:,"Revenue adjustment coefficient"])
    intensity_variable_cost = float.(land_use_inputs[:,"Variable expenses (\$/ha)"])
    nl_mitigation_variable_cost =  float.(land_use_inputs[:,:NL_Mitigation_Cost])
    methane_emission_rate = land_use_inputs[:,:Methane_emissions_factor]
    max_production_capability = land_use_inputs[:,"Production capability scalar"]

    #   Limit values (0% -> 100%) for intensity and NL_mitigation, set for each land_use type
    intensity_limits = [(0.0,1.0) for i in 1:length(land_use_names)]
    nl_mitigation_limits = [(0.0,1.0) for i in 1:length(land_use_names)]
    
    nl_mitigation_subsidy = 0.0

    #   Returns carbon sequestration values, if land_use supports
    carbon_sequestration_over_time = read_carbon_tables.(land_use_names)
    #   All zeroes on initialisation
    nitrous_oxide_price = 0.0
    methane_price = 0.0
    carbon_price = 0.0
    
    max_nl_mitigation_capability = 4.0
    production_capability_transform = float.(land_use_inputs[:,"Production capability transform"])
    
    #   Cost_switching for land_uses
    cost_to_switch = DataFrame(XLSX.readtable("data/Model inputs.xlsx","switching_cost_matrix"))
    cost_to_switch_from = read_cost_to_switch_matrix.(land_use_names,Ref(cost_to_switch))
    
    #   Spillover_learning for land_uses
    spillover_learning = DataFrame(XLSX.readtable("data/Model inputs.xlsx","transferability_matrix"))
    spillover_learning_rate_from = read_spillover_learning_matrix.(land_use_names,Ref(spillover_learning))

    return LandUse.(land_use_names,
                    land_use_categories,
                    fixed_costs,
                    years_to_yield,
                    crop_lifetime,
                    revenue_adjustment,
                    intensity_variable_cost,
                    nl_mitigation_variable_cost,
                    intensity_limits,
                    nl_mitigation_limits,
                    nl_mitigation_subsidy,
                    methane_emission_rate,
                    carbon_sequestration_over_time,
                    nitrous_oxide_price,
                    methane_price,
                    carbon_price,
                    max_production_capability,
                    max_nl_mitigation_capability,
                    production_capability_transform,
                    cost_to_switch_from,
                    spillover_learning_rate_from)
end

function test_profitability(row,land_use::LandUse)
    if land_use.name == "fallow" || occursin("forestry",land_use.name)
        return true
    end
    # average revenue over area of parcel
    revenue = land_use.revenue_adjustment * row["economic_indicator_revenue-"*land_use.name*"_value"] / (row[:total_area] / 10000)
    var_cost = land_use.intensity_variable_cost
    # test profit at high capability and optimal behaviour
    capability = 1.5 * land_use.max_production_capability
    opt_intensity = clamp(log(revenue * capability / var_cost) / capability,0,1)
    profit_at_optimum = revenue * (1 - exp(-capability * opt_intensity)) - (var_cost * opt_intensity) - land_use.fixed_costs
    return profit_at_optimum >= 0
end

function restrict_land_uses(land_uses,land_use_inputs,row;var_cost_threshold=0.0)
    allowed_land_uses::Vector{String} = []
    #   Iterate through each land use, check if revenue is either 0 or greather than or equal to total fixed costs (economic viability)
    #   Check if the average_angle is less than 5 degrees, if it is, exclude land-dependent uses
    non_irrigated_typecodes = [2,8,13,14,15,16,17,18]
    for lu in land_uses
        # disable dairy-irrigated for parcels with a non irrigable dairy typology
        if (lu.name == "dairy-irrigated") && (row["dairy_typology_typecode"] in non_irrigated_typecodes)
            continue
        end
        if "economic_indicator_revenue-"*lu.name*"_value" in names(row)
            if test_profitability(row,lu)
                max_slope = land_use_inputs[land_use_inputs[:,:Crop] .== lu.name,:Max_slope][1]
                if row["angl_avg"] <= max_slope
                    max_luc = string(land_use_inputs[land_use_inputs[:,:Crop] .== lu.name,:Max_LUC][1])
                    lu_cat = land_use_inputs[land_use_inputs[:,:Crop] .== lu.name,:Category][1]
                    existing_use = row[:seed_land_use]
                    # handle non-integer values in lriluc
                    lri = string(row["lri_luc"][1]) in [string(i) for i in 1:8] ? string(row["lri_luc"][1]) : "8"
                    # println(lri)
                    if (lri <= max_luc) || (lu_cat == existing_use)
                        push!(allowed_land_uses,lu.name)
                    end
                end
            end
        end
    end
    if length(allowed_land_uses) < 5
        println("Warning: ",row[:uid]," has ",length(allowed_land_uses)," allowed land uses")
        println("LUC: ",row["lri_luc"][1]," slope: ",row["angl_avg"]," degrees")
    end
    #   Returns allowed land uses - names of all land uses that are viable
    return allowed_land_uses
end
function create_land_parcels(land_uses,land_info)

    #   Initalise empty array to store land parcel objects
    parcels = []
    land_use_inputs = DataFrame(XLSX.readtable("data/Model inputs.xlsx","land_use_parameters"))
    dairy_typology_inputs = DataFrame(XLSX.readtable("data/Model inputs.xlsx","dairy_typology"))
    #   Iterate over each row of land_info, using enumerate to track row and index
    for (idx,row) in enumerate(eachrow(land_info))
        uid = row[:uid]
        nzseg = row[:catchment_nzsegment]
        #   Returns allowed land uses - names of all land uses that are viable
        lucas_2016_land_use = row[:seed_land_use]
        allowed_land_uses = restrict_land_uses(land_uses,land_use_inputs,row)

        land_use = "fallow"
        area_ha = row[:total_area]/10000
        avg_angle = row["angl_avg"] / area_ha
        years_in_current_land_use = 0.0
        
        #   Get land use parameters
        typecode, min_intensity_by_land_use, max_intensity_by_land_use, max_revenue_by_land_use, n_mitigation_potential_by_land_use, p_mitigation_potential_by_land_use, nitrogen_loss_coefficients_by_land_use, 
        phosphorous_loss_coefficients_by_land_use, nitrous_oxide_loss_coefficients_by_land_use, sediment_intensity_coefficients_by_land_use, sediment_loss_potenial_by_land_use, carbon_sequestration_adjustment_by_land_use = get_land_attributes(area_ha,allowed_land_uses,land_use_inputs,dairy_typology_inputs,row)
        
        
         
        nitrogen_price = 0.0
        phosphorous_price = 0.0
        sediment_price = 0.0
        intensity = 0.5
        nl_mitigation = 0.5
        utility = 0.0
        profit_npv_per_ha = 0.0
        total_yearly_profit = 0.0
        switching_costs_incurred = 0.0
        nl_mitigation_subsidy = 0.0
        wq_pollutant_costs = 0.0
        ghg_emission_costs = 0.0
        sequestration_payments = 0.0
        nitrogen_loss = 0.0
        phosphorous_loss = 0.0
        methane_emissions = 0.0
        nitrous_oxide_emissions = 0.0
        sediment_loss = 0.0
        carbon_sequestered = 0.0
        
        #   Create landparcel object with allowed land uses, parameters, and indicators for each land_info type
        parcel = LandParcel(uid,
                            nzseg,
                            lucas_2016_land_use,
                            allowed_land_uses, 
                            land_use, 
                            area_ha,
                            avg_angle,
                            years_in_current_land_use, 
                            max_revenue_by_land_use,  
                            typecode,
                            min_intensity_by_land_use,
                            max_intensity_by_land_use,
                            n_mitigation_potential_by_land_use,
                            p_mitigation_potential_by_land_use,
                            nitrogen_loss_coefficients_by_land_use,
                            phosphorous_loss_coefficients_by_land_use,
                            nitrous_oxide_loss_coefficients_by_land_use,
                            sediment_intensity_coefficients_by_land_use,
                            sediment_loss_potenial_by_land_use,
                            carbon_sequestration_adjustment_by_land_use,
                            nitrogen_price,
                            phosphorous_price,
                            sediment_price,
                            intensity,
                            nl_mitigation,
                            utility,
                            profit_npv_per_ha,
                            total_yearly_profit,
                            switching_costs_incurred,
                            nl_mitigation_subsidy,
                            wq_pollutant_costs,
                            ghg_emission_costs,
                            sequestration_payments,
                            nitrogen_loss,
                            phosphorous_loss,
                            methane_emissions,
                            nitrous_oxide_emissions,
                            sediment_loss,
                            carbon_sequestered
                            )

        #   Add landparcel object to list of objects for different land_info types
        push!(parcels,parcel)
    end
    return parcels
end

function merge_polygons(polygons)

    #   Creates empty multi-polygon
    multi_parcel_polygon = ArchGDAL.createmultipolygon()

    #   Iterates over each unit of multipolygon, checks it's type
    for p in polygons
        #   If single polygon, add directly to multipolygon container
         if typeof(p) == ArchGDAL.IGeometry{ArchGDAL.wkbPolygon}
            ArchGDAL.addgeom!(multi_parcel_polygon,p)
        #   Otherwise if a multi-polygon, iterate through individual polygons, adding each single polygon to container
         elseif typeof(p) == ArchGDAL.IGeometry{ArchGDAL.wkbMultiPolygon}
            for g in 0:ArchGDAL.getcoorddim(p)-1
                geom = ArchGDAL.getgeom(p,g)
                if typeof(geom) != ArchGDAL.IGeometry{ArchGDAL.wkbUnknown}
                    ArchGDAL.addgeom!(multi_parcel_polygon,geom)
                end
                
            end
        end
    end
#   Returns all merged single polygons
return multi_parcel_polygon
end


function get_agent_coords(land_parcel_ids,land_info)
    #   Returns all columns and rows of geometric for matched farms in land_info and land_parcel_ids
    polygons = filter(x->x.uid in land_parcel_ids,land_info)[:,:geom]

    #   Merges all polygons into a single multi-polygon
    multi_parcel_polygon = merge_polygons(polygons)

    #   Returns x and y coordinates of centroid of multi_parcel_polygon
    x_coordinate = ArchGDAL.getx(centroid(multi_parcel_polygon),0)
    y_coordinate = ArchGDAL.gety(centroid(multi_parcel_polygon),0)
    return x_coordinate,y_coordinate
end


function seed_initial_capability(land_uses,land_parcels,capability_initialisation_matrix)
    lucas_lus = unique([lp.lucas_2016_land_use for lp in land_parcels])
    p_a_capability_dict = Dict(lu.name=>0.0 for lu in land_uses)
    people_capability = rand(TruncatedNormal(0.4,0.25,0.001,1))
    business_capability = rand(TruncatedNormal(0.4,0.25,0.001,1))
    for lu in land_uses
        row = filter(x -> x[:land_use_category] == lu.category, capability_initialisation_matrix)
        capability_mean = maximum(eachcol(row[:,lucas_lus]))[1]
        p_a_capability_dict[lu.name] = rand(TruncatedNormal(capability_mean,0.15,0.001,1))
    end
    if "horticulture" in lucas_lus
        people_capability += 0.5
        business_capability += 0.5
    end
    return p_a_capability_dict, people_capability, business_capability
end


function seed_purchaser_capability(land_uses,land_parcels)
    capability_dict = Dict(lu.name=>0.0 for lu in land_uses)
    people_capability = rand(TruncatedNormal(0.4,0.25,0.001,1))
    business_capability = rand(TruncatedNormal(0.4,0.25,0.001,1))
    for lp in land_parcels
        init_cap = lp.lucas_2016_land_use == "horticulture" ? 0.9 : 0.6
        capability_dict[lp.land_use_name] = init_cap
    end
    for lu in land_uses
        capability_dict[lu.name] += rand(TruncatedNormal(0.2,0.15,0,0.5))
    end
    if any([lp.lucas_2016_land_use == "horticulture" for lp in land_parcels])
        people_capability += 0.5
        business_capability += 0.5
    end
    return capability_dict, people_capability, business_capability
end

function create_utility_barrier_matrix(inputs)
    #   Create empty dictionary to store utility barrier matrix
    # a dictionary of dictionaries is isomorphic to a matrix don't @ me
    utility_barrier_dict = Dict{String,Dict{String,Float64}}()
    for name in unique(inputs[:,:from])
        values = filter(row -> row[:from] == name, inputs)
        #   map column name from spill_over to first row in values, and store dictionary of values
        inner_dict = Dict{String,Float64}()
        for (k,v) in zip(names(inputs)[2:end],values[1,2:end])
            if (v == 1.0) || (v == 0.0)
                inner_dict[k] = v
            else
                # apply noise to positive utility values. Scale variance with sqrt value
                variance =  sqrt(v - 1) / 6
                inner_dict[k] = clamp(rand(LogNormal(v,variance)),0.8,4)
            end
        end
        utility_barrier_dict[name] = inner_dict
    end
    return utility_barrier_dict
end

function create_agent(uid, land_parcels, land_uses, x_coordinate, y_coordinate,capability_initialisation_matrix,utility_barrier_inputs)
    
    age = float(rand(20:64))
    # println(uid, " ", agent_land_parcels)
    likelihood_to_change = rand(TruncatedNormal(0.05,0.05,0.01,0.2)) * 0.99^(age-20)
    learning_rate = rand(TruncatedNormal(1,0.3,0.5,2))
    time_spent_on_external_resources = 0.0
    plant_and_animal_capability_by_land_use, people_capability, business_capability = capability_initialisation_matrix == "purchaser" ? seed_purchaser_capability(land_uses,land_parcels) : seed_initial_capability(land_uses,land_parcels,capability_initialisation_matrix)
    nl_mitigation_capability = rand(TruncatedNormal(0.1,0.3,0.001,1))
    environmental_concern = rand() < 0.3 ? 0 : rand(Uniform(0,10))
    nitrogen_weighting = environmental_concern * 3e-4
    phosphorous_weighting = environmental_concern * 5e-3
    ghg_weighting = environmental_concern * 1e-6
    sediment_weighting = environmental_concern * 1e-3
    financial_discount_rate = rand(Uniform(0.03,0.09))
    env_discount_rate = 0.5 * financial_discount_rate
    carbon_price_growth_rate_forecast = rand(TruncatedNormal(0.0,0.02,-0.03,0.04))
    deviation_parameter = rand(LogNormal(-0.1,0.1))
    utility_barrier_matrix = create_utility_barrier_matrix(utility_barrier_inputs)
    is_non_compliant = environmental_concern == 0 ? rand() < 0.67 : false
    profit_to_date = 0
    undertaking_initial_evaluation = true
    neighbours = []
    #   Create agent object each agent has a unique id, parameters, an initally empty list of neighbours, and attached land parcels associated 
    agent = Agent(uid,
                    age,
                    x_coordinate,
                    y_coordinate, 
                    likelihood_to_change, 
                    learning_rate, 
                    time_spent_on_external_resources,
                    plant_and_animal_capability_by_land_use,
                    people_capability,
                    business_capability, 
                    nl_mitigation_capability, 
                    nitrogen_weighting,
                    phosphorous_weighting,
                    ghg_weighting, 
                    sediment_weighting,
                    financial_discount_rate, 
                    env_discount_rate, 
                    carbon_price_growth_rate_forecast,
                    utility_barrier_matrix,
                    deviation_parameter,
                    is_non_compliant,
                    profit_to_date,
                    undertaking_initial_evaluation,
                    land_parcels,
                    neighbours)
    return agent
end


function create_agents(land_info, land_parcels, land_uses)
    #   Empty list for different agents
    agents = []
    capability_initialisation_matrix = DataFrame(XLSX.readtable("data/Model inputs.xlsx","capability_initialisation"))
    utility_barrier_inputs = DataFrame(XLSX.readtable("data/Model inputs.xlsx","utility_barrier_matrix"))
    #   For each unique farm ID
    for id in unique(land_info[:,:farm_id])
        uid = string("agent_", id)
        #   Create list of unique identifiers for land parcels that belong to specific farm 
        land_parcel_ids = [row[:uid] for row in eachrow(land_info) if row[:farm_id]==id]
        #   Returns x, y coordinates for centroid of land parcels
        x,y = get_agent_coords(land_parcel_ids,land_info)
        
        #   Select only those land parcels whose UIDs are in the land_parcel_ids list
        agent_land_parcels = filter(x->x.uid in land_parcel_ids,land_parcels)
        #   If no land parcels, skip
        if length(agent_land_parcels) == 0
            continue
        end
        #   Create agent object each agent has a unique id, parameters, an initally empty list of neighbours, and attached land parcels associated
        agent = create_agent(uid,agent_land_parcels,land_uses,x,y,capability_initialisation_matrix,utility_barrier_inputs)
        #   Add agent to list of agents
        push!(agents,agent)
    end
    return agents
end


function assign_slope_class(slope)
    if (slope === missing) | (slope < 5)
        return "Low"
    elseif slope < 15
        return "Moderate"
    else
        return "High"
    end
end

function get_existing_land_use(row)
    # Returns current land use based on LUCAS 2016 classification
    lucid = row["lucas_LUCID_2016"]
    if lucid === missing
        return "other"
    end
    if (lucid == 75) || (lucid == 76)
        subid = row["lucas_SUBID_2016"]
        if subid === missing
            return "other"
        end
        if subid == 502
            return "dairy"
        elseif subid == 503
            return "sheep-and-beef"
        else
            return "other"
        end
    elseif lucid == 77
        return "horticulture"
    elseif lucid == 78
        return "arable"
    else
        return "other"
    end
end

function read_and_process_parcels(geodata,sediment_model)
    # read in geospatial data
    gdf = GeoDataFrames.read(geodata)
    if sediment_model == "alt"
        gdf[:,"sediment_load_Herbaceous"] = gdf[:,"sediment_load_Herbaceous_alt"]
        gdf[:,"sediment_load_Trees_scrub"] = gdf[:,"sediment_load_Trees_scrub_alt"]
        gdf[:,"sediment_catchment_load_total_not_in_model"] = gdf[:,"sediment_catchment_load_total_not_in_model_alt"]
        gdf[:,"min_out_of_model_load"] = gdf[:,"min_out_of_model_load_alt"]
    end

    #Economic indicator columns
    gdf[:,"economic_indicator_revenue-fallow_value"] = zeros(nrow(gdf))
    gdf[:,"economic_indicator_revenue-carbon-forestry_value"] = zeros(nrow(gdf))
    gdf[:,"economic_indicator_revenue-native-forestry_value"] = zeros(nrow(gdf))
    gdf[:,"economic_indicator_revenue-production-forestry_value"] = gdf[:,"economic_indicator_avg_Radiata300_m3_per_ha"] .* 45000 / 27
    
    # FIX with updated pasture yield values from Giotto
    # gdf[:,"economic_indicator_revenue-dairy_value"] = gdf[:,"DairyNZ revenue per ha"]
    # gdf[:,"economic_indicator_revenue-dairy-irrigated_value"] = gdf[:,"DairyNZ revenue per ha"]
    # gdf[:,"economic_indicator_revenue-sheep-and-beef_value"] = gdf[:,"B+LNZ revenue per ha"]

    # convert pasture yield to median revenue equivalent - mean rev / mean pasture yield from calibration against previous initialisation run
    dairy_pasture_yield_scalar = 9454/17724
    snb_pasture_yield_scalar = 894/12431
    gdf[:,"economic_indicator_revenue-dairy_value"] = gdf[:,"PastureYield_Rainfed200NCap_25m"] * dairy_pasture_yield_scalar
    gdf[:,"economic_indicator_revenue-dairy-irrigated_value"] = gdf[:,"PastureYield_Irrigated200NCap_25m"] * dairy_pasture_yield_scalar
    gdf[:,"economic_indicator_revenue-sheep-and-beef_value"] = gdf[:,"PastureYield_Rainfed200NCap_25m"] * snb_pasture_yield_scalar
    
    
    
    gdf[:,"total_area"] = 10000 .* ones(nrow(gdf))
    gdf[:,:seed_land_use] = get_existing_land_use.(eachrow(gdf))
    # assign_slope_class function applied for slope class column
    gdf[:,:slope_class] = map(assign_slope_class,gdf[:,"angl_avg"])

    #  set of catagorical value columns
    idcols = ["lri_luc","farm_id","seed_land_use","slope_class","dairy_typology_typecode","catchment_nzsegment","catchment_REC2_TerminalSegment","catchment_WaterShedAreaKM2","catchment_MAL_Nitrogen","catchment_MAL_Phosphorous","catchment_MAL_Sediment","sediment_catchment_load_total_not_in_model","min_out_of_model_load"]

    #  set of numerical value columns
    numcols = ["total_area","sediment_load_Herbaceous","sediment_load_Trees_scrub","angl_avg"]

    #  Iterate over columns of gdf dataframe, check if column is revenue indicator, and retrieve data frame of revenue indicators
    eco_indicators = [c for c in names(gdf) if occursin("economic_indicator_revenue",c)]

    #  Missing numeric values set to 0
    #  Columns that are fully empty are removed
    for e in vcat(eco_indicators,numcols,idcols)
        if eltype(gdf[:,e]) == Union{Missing, Float64}
            gdf[:,e] = replace(gdf[:,e],missing=>0.0)
        elseif eltype(gdf[:,e]) == Union{Missing, String}
            gdf[:,e] = replace(gdf[:,e],missing=>"NA")
        elseif eltype(gdf[:,e]) == Missing
            # handle fully-missing column (really shouldn't happen but does)
            select!(gdf, Not(e))
            gdf[:,e] = [0.0 for i in 1:nrow(gdf)]
        # else
        #     println("Warning: ",e," is of type ",eltype(gdf[:,e]))
        #     # gdf[:,e] = [0.0 for i in 1:nrow(gdf)]
        end
    end

    #  Remove rows where a value in any column is missing
    newdf = disallowmissing!(gdf[:,vcat(idcols,numcols,["geom"],eco_indicators)])

    #   Farms are in individual hectares data form -> Group by same overall farm (based on [idcols]) break into subfarms based on sediment and slop differences (idcols)
    merged = combine(groupby(newdf,idcols),:geom=>merge_polygons, eco_indicators .=> sum, numcols .=> sum, renamecols=false)
    # re-average angle instead of summing.
    merged[:,"angl_avg"] = merged[:,"angl_avg"] ./ (merged[:,"total_area"]/10000)
    return merged
end


function create_catchments(land_parcels,land_info)
    catchments::Array{Catchment,1} = []
    for lp in land_parcels
        if lp.nzseg in [c.nzseg for c in catchments]
            # add land parcel to existing catchment
            catchment = catchments[findfirst(c->c.nzseg==lp.nzseg,catchments)]
            push!(catchment.land_parcels,lp)
            catchment.area_in_forestry_or_fallow += lp.area_ha

        else
            # create new catchment
            nzseg = lp.nzseg
            land_parcels = [lp]
            primary_lu = "fallow"
            tseg = land_info[land_info[:,:catchment_nzsegment] .== lp.nzseg,:catchment_REC2_TerminalSegment][1]
            watershed = land_info[land_info[:,:catchment_nzsegment] .== lp.nzseg,:catchment_WaterShedAreaKM2][1] * 100
            area_productive = 0.0
            area_forestry_or_fallow = lp.area_ha
            total_coverage = 0.0
            MAL_N = land_info[land_info[:,:catchment_nzsegment] .== lp.nzseg,:catchment_MAL_Nitrogen][1]
            MAL_P = land_info[land_info[:,:catchment_nzsegment] .== lp.nzseg,:catchment_MAL_Phosphorous][1]
            MAL_Sediment = land_info[land_info[:,:catchment_nzsegment] .== lp.nzseg,:catchment_MAL_Sediment][1]
            out_of_model_Sediment_load = land_info[land_info[:,:catchment_nzsegment] .== lp.nzseg,:sediment_catchment_load_total_not_in_model][1]
            min_out_of_model_Sediment_load = land_info[land_info[:,:catchment_nzsegment] .== lp.nzseg,:min_out_of_model_load][1]
            has_Sediment_exemption = false
            Load_N = 0.0
            Load_P = 0.0
            Load_Sediment = 0.0
            excess_N_percent = 0.0
            excess_P_percent = 0.0
            excess_Sediment_percent = 0.0
            tax_revenue_pool = 0.0
            
            catchment = Catchment(nzseg,
                                tseg,
                                land_parcels,
                                primary_lu,
                                watershed,
                                area_productive,
                                area_forestry_or_fallow,
                                total_coverage,
                                MAL_N,
                                MAL_P,
                                MAL_Sediment,
                                out_of_model_Sediment_load,
                                min_out_of_model_Sediment_load,
                                has_Sediment_exemption,
                                Load_N,
                                Load_P,
                                Load_Sediment,
                                excess_N_percent,
                                excess_P_percent,
                                excess_Sediment_percent,
                                tax_revenue_pool)
            push!(catchments,catchment)
        end
    end
    return catchments
end


function change_params_for_SA!(sensitivity_analysis,agents,land_uses)
    paramtable = DataFrame(XLSX.readtable("calibration/sensitivity analysis.xlsx","Sheet1"))
    row = filter(row -> row[:sensitivity_analysis_scenario_number] == sensitivity_analysis, paramtable)[1,:]
    if !ismissing(row[:likelihood_to_change])
        for agent in agents
            agent.likelihood_to_change = row[:likelihood_to_change]
        end
    end
    if !ismissing(row[:learning_rate])
        for agent in agents
            agent.learning_rate = row[:learning_rate]
        end
    end
    if !ismissing(row[:capability])
        for agent in agents
            agent.people_capability = row[:capability]
            agent.business_capability = row[:capability]
            agent.nl_mitigation_capability = row[:capability]
            for k in keys(agent.plant_and_animal_capability_by_land_use)
                agent.plant_and_animal_capability_by_land_use[k] = row[:capability]
            end
        end
    end
    if !ismissing(row[:environmental_concern])
        for agent in agents
            agent.nitrogen_weighting = row[:environmental_concern] * 3e-4
            agent.phosphorous_weighting = row[:environmental_concern] * 5e-3
            agent.ghg_weighting = row[:environmental_concern] * 1e-6
            agent.sediment_weighting = row[:environmental_concern] * 1e-3
        end
    end
    if !ismissing(row[:deviation_parameter])
        for agent in agents
            agent.deviation_parameter = row[:deviation_parameter]
        end
    end
    if !ismissing(row[:financial_discount_rate])
        for agent in agents
            agent.financial_discount_rate = row[:financial_discount_rate]
            agent.env_discount_rate = 0.5 * row[:financial_discount_rate]
        end
    end
    if !ismissing(row[:carbon_price_growth_rate_forecast])
        for agent in agents
            agent.carbon_price_growth_rate_forecast = row[:carbon_price_growth_rate_forecast]
        end
    end
    if !ismissing(row[:utility_barrier])
        for agent in agents
            for k in keys(agent.utility_barrier_matrix)
                for k2 in keys(agent.utility_barrier_matrix[k])
                    agent.utility_barrier_matrix[k][k2] = row[:utility_barrier]
                end
            end
        end
    end
    if !ismissing(row[:nl_mitigation_cost])
        for lu in land_uses
            lu.nl_mitigation_variable_cost = row[:nl_mitigation_cost]
            if lu.category == "sheep-and-beef"
                lu.nl_mitigation_variable_cost = row[:nl_mitigation_cost] * 0.2
            end
        end
    end
end

function generate_initial_population(geodata,max_parcels,sediment_model)
    
    #   Create land groups for different farms
    land_info = @time read_and_process_parcels(geodata,sediment_model)
    if nrow(land_info) > max_parcels
        land_info = land_info[1:max_parcels,:]
    end

    #   Label land_info rows in uid column, land_parcel_n
    land_info[:,:uid] = [string("land_parcel_", idx) for idx in 1:nrow(land_info)]

    #   Initialise the different land uses and their associated values
    land_uses = @time initialise_land_uses()

    #   Using the land_info & land_uses, create land parcel objects which connect land area to usage parameters
    land_parcels = @time create_land_parcels(land_uses,land_info)

    #   Variable agents should be an array of Agent objects, created with create_agents
    agents::Array{Agent,1} = @time create_agents(land_info, land_parcels, land_uses)

    catchments::Array{Catchment,1} = create_catchments(land_parcels,land_info)

    println(length(agents)," agents managing ",length(land_parcels)," parcels across ",length(catchments)," catchments")
    return agents, land_uses, land_info, catchments
end


function adjust_l2c!(agents)
    for agent in agents
        snb_area = sum([lp.area_ha for lp in agent.land_parcels if occursin("sheep-and-beef",lp.land_use_name)])
        total_area = sum([lp.area_ha for lp in agent.land_parcels])
        if snb_area > 0.5 * total_area
            agent.likelihood_to_change *= 0.1
        end
    end
end
end