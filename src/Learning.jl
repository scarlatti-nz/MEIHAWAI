module Learning
using ..Structs, Distributions, StatsBase, Random


function apply_spillover_learning!(agent::Agent,land_uses::Vector{LandUse},primary_land_use_name,capability_shift)
    #   Iterate over each land_use
    for lu in land_uses
        #   Updates the agent's plant and animal capability for each land use
        #   Adds (capability shift * spillover learning rate from the primary land use) to existing plant and animal capability for the land use
        agent.plant_and_animal_capability_by_land_use[lu.name] += capability_shift * lu.spillover_learning_rate_from[primary_land_use_name]
    end
end

function get_lu_weights(agent::Agent,neighbour::Agent,land_uses::Vector{LandUse})
    land_use_weights::Vector{Float64} = []
    for lu in land_uses
        if any([lp.land_use_name == lu.name for lp in intersect(agent.land_parcels,neighbour.land_parcels)])
            push!(land_use_weights,1)
        elseif any([lp.land_use_name == lu.name for lp in union(agent.land_parcels,neighbour.land_parcels)])
            push!(land_use_weights,0.5)
        else
            push!(land_use_weights,0.1)
        end
    end
    return land_use_weights
end

function conversation!(agent,neighbour,land_uses,primary_land_use_name,topic,capability_shift)
    fieldname = Symbol(topic*"_capability_by_land_use")
    agent_cap_dict = getfield(agent,fieldname)
    neighbour_cap_dict = getfield(neighbour,fieldname)
    if agent_cap_dict[primary_land_use_name] < neighbour_cap_dict[primary_land_use_name]
        capability_shift = (neighbour_cap_dict[primary_land_use_name] - agent_cap_dict[primary_land_use_name]) * capability_shift * agent.learning_rate
        apply_spillover_learning!(agent,land_uses,primary_land_use_name,capability_shift)
    else
        capability_shift = (agent_cap_dict[primary_land_use_name] - neighbour_cap_dict[primary_land_use_name]) * capability_shift * neighbour.learning_rate
        apply_spillover_learning!(neighbour,land_uses,primary_land_use_name,capability_shift)
    end
end

function conversation!(agent,neighbour,topic,capability_shift)
    fieldname = Symbol(topic*"_capability")
    agent_cap = getfield(agent,fieldname)
    neighbour_cap = getfield(neighbour,fieldname)
    if agent_cap < neighbour_cap
        agent_cap = agent_cap + (neighbour_cap - agent_cap) * capability_shift * agent.learning_rate
        setfield!(agent,fieldname,agent_cap)
    else
        neighbour_cap = neighbour_cap + (agent_cap - neighbour_cap) * capability_shift * neighbour.learning_rate
        setfield!(neighbour,fieldname,neighbour_cap)
    end
end

function learn_from_neighbours!(agent::Agent, land_uses::Vector{LandUse})
    # calibrated on workshop responses based on estimate of 100 conversations per year. doubled to account for skill discrepancy
    capability_shift = 0.0006
    expected_meetings_per_year = 10
    land_uses_not_fallow = [lu for lu in land_uses if lu.name != "fallow"]
    for neighbour in agent.neighbours
        number_of_meetings = rand(Poisson(expected_meetings_per_year))
        land_use_weights = get_lu_weights(agent,neighbour,land_uses_not_fallow)
        for meeting in 1:number_of_meetings
            topic = sample(["plant_and_animal","people","business","nl_mitigation"],1)[1]
            if topic == "plant_and_animal"
                land_use = sample(land_uses_not_fallow,Weights(land_use_weights),1)[1]
                conversation!(agent,neighbour,land_uses,land_use.name,topic,capability_shift)
            else
                conversation!(agent,neighbour,topic,capability_shift)
            end
        end
    end
end

function learn_by_doing!(agent,land_uses;learning_activity_threshold=0.01)
    #   Filters out land uses not "fallow"
    land_uses_not_fallow = [lu for lu in land_uses if lu.name != "fallow"]
   
    #   Creates list of the area of each land parcel owned by the agent if the land use of the parcel is not "fallow"
    #   Sums list to find total area, clamps, so value is set to lower or upper limit if unbound
    total_area_farmed = clamp(sum([lp.area_ha for lp in agent.land_parcels if lp.land_use_name in [lu.name for lu in land_uses_not_fallow]]),1,Inf)
    
    # values hard coded from workshop responses spreadsheet
    capability_shift_plant_and_animal = 0.0049
    capability_shift_people = 0.0021
    capability_shift_business = 0.0012
    capability_shift_nl_mitigation = 0.0021

    #   Updates agents people and business capability based on hours, learning rate, and quality
    agent.people_capability +=  agent.learning_rate * capability_shift_people
    agent.business_capability +=  agent.learning_rate * capability_shift_business

    #   Check if agent's land parcels have non-linear mitigation greater than the learning activity threshold
    #   Boolean list - true: nl_mitigation > learning_threshold, false: otherwise
    if any([lp.nl_mitigation > learning_activity_threshold for lp in agent.land_parcels])
        #   Updates agent's nl_mitigation_capability based on hours, learning rate, and quality if boolean is true
        agent.nl_mitigation_capability +=  agent.learning_rate * capability_shift_nl_mitigation
    end

    #   Iterate over each land use that is not "fallow"
    for lu in land_uses_not_fallow
        #   Calculates the hours spent on the current land use
        hours_spent = sum([lp.area_ha for lp in agent.land_parcels if lp.land_use_name == lu.name]) / total_area_farmed
        #   Calculate the plant and animal shift = Improvemenet of agents understand of plant and animal capability 
        plant_and_animal_shift = hours_spent * agent.learning_rate * capability_shift_plant_and_animal
        
        #   agent learns by doing each land use and applies spillover learning to the agent
        apply_spillover_learning!(agent,land_uses,lu.name,plant_and_animal_shift)
    end
end

function get_participant_weighting(agents,resource)
    # agents cannot spend more than 100 hours per year engaging with external resources
    resource_time_cap = 100
    l2c_weighting = [a.likelihood_to_change for a in agents]
    primary_topic = sort(collect(resource.capability_shift_by_topic),by=x->x[2],rev=true)[1][1]
    saturation_weighting = [a.time_spent_on_external_resources + resource.duration > resource_time_cap ? 0 : 1 for a in agents]
    if primary_topic == "plant_and_animal"
        capability_weighting = [max(0,1 - (0.8 - maximum(values(a.plant_and_animal_capability_by_land_use)))^2) for a in agents]
    else
        # quadratic weighting, most likely to engage are those roughly in the middle of the spectrum
        capability_weighting = [max(0,1 - (0.8 - getfield(a,Symbol(primary_topic*"_capability")))^2) for a in agents]
    end
    return Weights(l2c_weighting .* capability_weighting .* saturation_weighting)
end

function learn_from_resources!(agents::Vector{Agent},land_uses::Vector{LandUse},external_resources::Vector{ExternalResource})
    # randomise the order so that when hitting the cap, the same resource is not always discarded
    for resource in shuffle(external_resources)
        weights = get_participant_weighting(agents,resource)

        n_participants = min(resource.n_participants,length([w for w in weights if w>0]))
        if n_participants == 0
            continue
        end
        resource.run_this_timestep = true
        participants = sample(agents,weights,n_participants,replace=false)
        for participant in participants
            participant.time_spent_on_external_resources += resource.duration
            for topic in ["people","business","nl_mitigation"]
                participant_capability = getfield(participant,Symbol(topic*"_capability"))
                participant_capability += resource.capability_shift_by_topic[topic] * participant.learning_rate
                setfield!(participant,Symbol(topic*"_capability"),participant_capability)
            end
            # get the land use for which the agent has the greatest actively farmed area
            minus_fallow_and_forestry = [lu for lu in land_uses if lu.name != "fallow" && !occursin("forestry",lu.name)]
            lu_by_area = Dict(lu.name => sum([lp.area_ha for lp in participant.land_parcels if lp.land_use_name == lu.name]) for lu in minus_fallow_and_forestry)
            if sum(values(lu_by_area)) == 0
                continue
            end
            primary_land_use = sort(collect(lu_by_area),by=x->x[2],rev=true)[1][1]
            if rand() > resource.land_use_change_focus
                capability_shift = resource.capability_shift_by_topic["plant_and_animal"] * participant.learning_rate
                apply_spillover_learning!(participant,land_uses,primary_land_use,capability_shift)
            else
                arable_hort_only = [lu for lu in minus_fallow_and_forestry if !occursin("dairy",lu.name) && !occursin("sheep-and-beef",lu.name)]
                if length(arable_hort_only) == 0
                    continue
                end
                alt_land_use = sample([lu.name for lu in arable_hort_only if lu.name != primary_land_use],1)[1]
                capability_shift = resource.capability_shift_by_topic["plant_and_animal"] * participant.learning_rate
                apply_spillover_learning!(participant,land_uses,alt_land_use,capability_shift)
            end
        end
    end
end

end