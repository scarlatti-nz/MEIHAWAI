module Policies
using ..Structs, ..Population, StatsBase, Distributions, DataFrames, XLSX

function increase_wq_prices!(agents::Vector{Agent})
        # updates local nutrient prices with 3% growth rate
    growth_rate = 1.03
    for a in agents
        for lp in a.land_parcels
            for fieldname in fieldnames(LandParcel)
                field = getfield(lp,fieldname)
                if occursin("price",string(fieldname))
                    field *= growth_rate
                    setfield!(lp,fieldname,field)
                end
            end
        end
    end
end

function exponential_price_path!(land_uses::Vector{LandUse})
    growth_rate = 1.03
    # updates ghg emissions and sequestration prices with 3% growth rate
    for lu in land_uses
        for fieldname in fieldnames(LandUse)
            field = getfield(lu,fieldname)
            if occursin("price",string(fieldname))
                field *= growth_rate
                setfield!(lu,fieldname,field)
            end
        end
    end
end

function parabolic_price_path!(land_uses::Vector{LandUse},timestep,N2OPrice,CH4Price,CO2Price)
    # updates ghg emissions and sequestration prices with 3% growth rate
    quadratic_term = -0.000689
    linear_term = 0.0342
    # remember to divide by 1000 because prices are in $/tonne and we want $/kg
    price_coefficient_at_time_t = (1 + linear_term * timestep + quadratic_term * timestep^2) / 1000
    for lu in land_uses
        if occursin("forestry",lu.name)
            lu.carbon_price = price_coefficient_at_time_t * CO2Price
        else
            lu.methane_price = price_coefficient_at_time_t * CH4Price
            lu.nitrous_oxide_price = price_coefficient_at_time_t * N2OPrice
        end
    end
end

function enable_bolus!(land_uses::Vector{LandUse},agents::Vector{Agent},enable_bolus)
    if enable_bolus
        #   Adds additional land_uses to vector if bolus enabled
        new_lus = Population.initialise_land_uses(enable_additional_uses=["dairy-bolus","sheep-and-beef-bolus"])
        for lu in new_lus
            push!(land_uses,lu)
        end

        #   Iterate over each agent and each field of agent, if dict with string keys and float64 values, add new land_uses assign same values as dairy and sheep-and-beef
        for agent in agents
            for fieldname in fieldnames(Agent)
                field = getfield(agent,fieldname)
                if field isa Dict{String,Float64}
                    field["dairy-bolus"] = field["dairy"]
                    field["sheep-and-beef-bolus"] = field["sheep-and-beef"]
                end
            end

            #   Iterate over land_parcels and add bolus uses to allowed_land_uses. Iterate over each field in land_parcel type, assign the same values as dairy and sheep-and-beef 
            for lp in agent.land_parcels
                if "dairy" in lp.allowed_land_uses
                    push!(lp.allowed_land_uses,"dairy-bolus")
                    for fieldname in fieldnames(LandParcel)
                        field = getfield(lp,fieldname)
                        if field isa Dict{String,Float64}
                            field["dairy-bolus"] = field["dairy"]
                        end
                    end
                end
                if "sheep-and-beef" in lp.allowed_land_uses
                    push!(lp.allowed_land_uses,"sheep-and-beef-bolus")
                    for fieldname in fieldnames(LandParcel)
                        field = getfield(lp,fieldname)
                        if field isa Dict{String,Float64}
                            field["sheep-and-beef-bolus"] = field["sheep-and-beef"]
                        end
                    end
                end
            end
        end
    end
end


function calculate_local_N_price(catchment::Catchment,NPrice;min_effective_price_per_ha = 100,max_effective_price_per_ha = 10000)
    if (NPrice == 0) || (catchment.Load_N <= 0) || (catchment.MAL_N <= 0) || (catchment.Load_N < catchment.MAL_N)
        return 0
    end
    total_N_liability =  NPrice * (catchment.Load_N - catchment.MAL_N)
    local_N_price_per_kg = total_N_liability / catchment.Load_N
    N_loss_worst_ha = maximum([lp.nitrogen_loss for lp in catchment.land_parcels])
    local_N_price_per_kg = clamp(local_N_price_per_kg, min_effective_price_per_ha / N_loss_worst_ha, max_effective_price_per_ha / N_loss_worst_ha)
    return local_N_price_per_kg
end

function calculate_local_P_price(catchment::Catchment,PPrice;min_effective_price_per_ha = 100,max_effective_price_per_ha = 10000)
    if (PPrice == 0) || (catchment.Load_P <= 0) || (catchment.MAL_P <= 0) || (catchment.Load_P < catchment.MAL_P)
        return 0
    end
    total_P_liability =  PPrice * (catchment.Load_P - catchment.MAL_P)
    local_P_price_per_kg = total_P_liability / catchment.Load_P
    P_loss_worst_ha = maximum([lp.phosphorous_loss for lp in catchment.land_parcels])
    local_P_price_per_kg = clamp(local_P_price_per_kg, min_effective_price_per_ha / P_loss_worst_ha, max_effective_price_per_ha / P_loss_worst_ha)
    return local_P_price_per_kg
end

function calculate_local_sediment_price(catchment::Catchment,SedimentPrice;min_effective_price_per_ha = 100,max_effective_price_per_ha = 10000)
    if (SedimentPrice == 0) || (catchment.Load_Sediment <= 0) || (catchment.MAL_Sediment <= 0) || (catchment.Load_Sediment < catchment.MAL_Sediment) || (catchment.has_Sediment_exemption)
        return 0
    end
    total_Sediment_liability =  SedimentPrice * (catchment.Load_Sediment - catchment.MAL_Sediment)
    local_Sediment_price_per_ton = total_Sediment_liability / catchment.Load_Sediment
    sediment_loss_worst_ha = maximum([lp.sediment_loss for lp in catchment.land_parcels])
    local_Sediment_price_per_ton = clamp(local_Sediment_price_per_ton, min_effective_price_per_ha / sediment_loss_worst_ha, max_effective_price_per_ha / sediment_loss_worst_ha)
    return local_Sediment_price_per_ton
end


function apply_incentives!(land_uses::Vector{LandUse},agents::Vector{Agent},catchments::Vector{Catchment},NPrice,PPrice,SedimentPrice,N2OPrice,CH4Price,CO2Price)
    #   Iterate over each land, use, if forestry, change nl_price and ghg_price
    for lu in land_uses
        if !occursin("forestry",lu.name)
            # convert from price per ton to price per kg here
            lu.methane_price = CH4Price / 1000
            lu.nitrous_oxide_price = N2OPrice / 1000
        else
            lu.carbon_price = CO2Price / 1000
        end
    end

    # don't set prices if MAL is 0, those are bad data cells.
    for c in catchments
        catchment_forest_sediment_production = sum([lp.sediment_loss_potenial_by_land_use["native-forestry"] * lp.area_ha for lp in c.land_parcels])
        if (c.min_out_of_model_Sediment_load + catchment_forest_sediment_production) > c.MAL_Sediment
            c.has_Sediment_exemption = true
        end
        local_N_price_per_kg = calculate_local_N_price(c,NPrice)
        local_P_price_per_kg = calculate_local_P_price(c,PPrice)
        local_Sediment_price_per_ton = calculate_local_sediment_price(c,SedimentPrice)

        for lp in c.land_parcels
            lp.nitrogen_price = local_N_price_per_kg
            lp.phosphorous_price = local_P_price_per_kg
            lp.sediment_price = local_Sediment_price_per_ton
        end
    end
end


### DEPRECATED
function apply_regulations!(land_uses::Vector{LandUse},agents::Vector{Agent},regulation_scale)

    #   Reads in data frame of regulations
    inputs = DataFrame(XLSX.readtable("data/Model inputs.xlsx","regulations"))

    #   Extracts default values from data frame for each regulation
    intensity_limit = inputs[inputs.Regulation .== "Restrictions on intensity", :default_value][1]
    nl_mitigation_limit = inputs[inputs.Regulation .== "Minimum mitigation requirement", :default_value][1]
    afforestation_percentile = inputs[inputs.Regulation .== "Mandatory afforestation of less productive land", :default_value][1]

    #   If regulation scale is greater than 0, iterate over each land use. Create array with intensities of land parcels of all agents with the same land_use
    if regulation_scale > 0
        for lu in land_uses
            intensity_distribution = [lp.intensity for a in agents for lp in a.land_parcels if lp.land_use_name == lu.name]
            if length(intensity_distribution)<10
                continue
            end
            #   Calculates the 90th percentile of intensity distribution
            upper_limit = quantile(intensity_distribution,(1 - (1 - intensity_limit/100) * regulation_scale))
            
            #   Sets bounds for intensity and nutrient loss mitigation
            lu.intensity_limits = (lu.intensity_limits[1], upper_limit)
            lu.nl_mitigation_limits = (nl_mitigation_limit * regulation_scale, lu.nl_mitigation_limits[2])
        end
        #   Sort parcels by sediment loss by land use in descending order and set the allowed land uses for the top n_afforestation parcels to forestry
        all_parcels = [lp for agent in agents for lp in agent.land_parcels]
        sort!(all_parcels,by=x->x.sediment_loss_potenial_by_land_use["dairy"],rev=true)
        n_parcels = length(all_parcels)
        n_afforestation = round(Integer,n_parcels * afforestation_percentile/100 * regulation_scale)
        for i in 1:n_afforestation
            all_parcels[i].allowed_land_uses = [lu.name for lu in land_uses if occursin("forestry",lu.name)]
        end
    end
end


function generate_baseline_external_resources(agents::Vector{Agent})::Vector{ExternalResource}
    external_resources = []
    topics = ["plant_and_animal","people","business","nl_mitigation"]

    #   Reads in dataframe from Excel file in extension sheet
    inputs = DataFrame(XLSX.readtable("data/Model inputs.xlsx","extension"))

    inputs = filter(row -> row[:Scenario] == "Baseline", inputs) 
    #   counts the number of agents
    n_agents = length(agents)

    # allocate extension spending to this grid cell based on national total of 200M and 68325 agents
    total_spend_available = 2e8 * n_agents / 68325

    #   Iterates over each row in inputs
    for row in eachrow(inputs)
        # don't want to read n_events_per timestep from the spreadsheet, we want to calculate this on the fly based on cost
        n_events = total_spend_available * row[:proportion_of_spend] / row[:cost_per_event]
        
        #   For each event:
        for i in 1:round(Integer,n_events)
            #   Assign event ID
            uid = string(row["Event type"], i)
            funding = "Baseline"
            capability_shift_by_topic = Dict(t=>row["capability_shift_"*t] for t in topics)
            n_participants = row[:n_participants_per_event]
            duration = row[:duration]
            land_use_change_focus = row[:land_use_change_focus]
            cost_per_event = row[:cost_per_event]
            run_this_timestep = false
            resource = ExternalResource(uid,funding,capability_shift_by_topic,n_participants,duration,land_use_change_focus,cost_per_event,run_this_timestep)
            push!(external_resources,resource)
        end
    end
    return external_resources
end

function generate_additional_external_resources(agents::Vector{Agent},extension_scale)::Vector{ExternalResource}
    external_resources = []

    #   List of topics 
    topics = ["plant_and_animal","people","business","nl_mitigation"]

    #   Reads in dataframe from Excel file in extension sheet
    inputs = DataFrame(XLSX.readtable("data/Model inputs.xlsx","extension"))

    inputs = filter(row -> row[:Scenario] == "Additional", inputs) 
    #   counts the number of agents
    n_agents = length(agents)

    # allocate extension spending to this grid cell based on national total of 200M and 68325 agents
    total_spend_available = extension_scale * 2e8 * n_agents / 68325

    #   Iterates over each row in inputs
    for row in eachrow(inputs)
        # don't want to read n_events_per timestep from the spreadsheet, we want to calculate this on the fly based on cost
        n_events = total_spend_available * row[:proportion_of_spend] / row[:cost_per_event]
        
        #   For each event:
        for i in 1:round(Integer,n_events)
            #   Assign event ID
            uid = string(row["Event type"], i)
            funding = "Additional"
            capability_shift_by_topic = Dict(t=>row["capability_shift_"*t] for t in topics)
            n_participants = row[:n_participants_per_event]
            duration = row[:duration]
            land_use_change_focus = row[:land_use_change_focus]
            cost_per_event = row[:cost_per_event]
            run_this_timestep = false
            resource = ExternalResource(uid,funding,capability_shift_by_topic,n_participants,duration,land_use_change_focus,cost_per_event,run_this_timestep)
            push!(external_resources,resource)
        end
    end
    return external_resources
end

end