module ModelStep

using ..Structs, ..Learning, ..Population, ..DecisionProblem, ..Policies, ..GenerateOutcomes, ..SocialNetwork, Distributions, DataFrames, XLSX

function create_successor_agent(agent,land_uses)
    # actually don't need to create a new agent, just reset age and likelihood to change and perturb capabilities 
    agent.age = rand(18:35)
    agent.likelihood_to_change = rand(TruncatedNormal(0.05,0.05,0.01,0.2)) * 0.99^(agent.age-20)
    for lu in land_uses
        old_cap = agent.plant_and_animal_capability_by_land_use[lu.name]
        agent.plant_and_animal_capability_by_land_use[lu.name] = rand(TruncatedNormal(old_cap,0.1,0.1,old_cap + 0.2))
    end
    agent.people_capability = rand(TruncatedNormal(agent.people_capability,0.1,0.1,agent.people_capability+0.2))
    agent.business_capability = rand(TruncatedNormal(agent.business_capability,0.1,0.1,agent.business_capability+0.2))
    agent.nl_mitigation_capability = rand(TruncatedNormal(agent.nl_mitigation_capability,0.1,0.1,agent.nl_mitigation_capability+0.2))
    agent.profit_to_date = 0.0
    return agent
end

function create_purchaser_agent(agent,land_uses,utility_barrier_inputs)
    new_agent = Population.create_agent(agent.uid,agent.land_parcels,land_uses,agent.x_coordinate,agent.y_coordinate,"purchaser",utility_barrier_inputs)
    # println("Purchaser agent created ",new_agent.uid)
    # purchaser agent behaves as if they are starting from scratch, acting as if coming to fallow land and evaluating all options equally.
    for lp in new_agent.land_parcels
        lp.years_in_current_land_use = 0.0
        lp.land_use_name == "fallow"
    end
    return new_agent
end

function retire(agent,land_uses,utility_barrier_inputs)

    #    Create a new agent using population.create_agent and returns agent object
    succession_rate = 0.5 # 50% (https://www.nzherald.co.nz/the-country/news/future-of-farming-communication-vital-in-farm-succession-planning/TB6IZ72NY3VNGKKAT5SDNMFJJQ/) chance of being replaced by a successor, who inherits capability and other values with small random perturbation and continues farm operations.
    inherit_neighbours = agent.neighbours
    if rand() < succession_rate
        new_agent = create_successor_agent(agent,land_uses)
    else
        new_agent = create_purchaser_agent(agent,land_uses,utility_barrier_inputs)
    end
    new_agent.neighbours = inherit_neighbours
    return new_agent
end

function add_initial_capability!(agent,land_uses)
    #   Set up intial capability for each land use for the time horizion, based on land_use before the time horizion 

    #   Calculate years each land parcel of agent has spent in current land use, excluding fallow
    years_spent_by_parcel = [agent.age - 18 for lp in agent.land_parcels if lp.land_use_name != "fallow"]

    #   If an agent has spent more than one year in a land use, they gain capability in that land use
    if !isempty(years_spent_by_parcel) && (maximum(years_spent_by_parcel) > 1)
        
        #   For land use, the agent's capability is increased each year spent in that land use
        for i in 1:maximum(years_spent_by_parcel)
            Learning.learn_by_doing!(agent,land_uses)
        end
    end
end

function advance_time(agent::Agent,land_uses::Vector{LandUse},timedelta,utility_barrier_inputs::DataFrame)
    #   Increment age and decrease likelihood to change land use
    agent.age += timedelta

    agent.likelihood_to_change *= 0.99
    agent.time_spent_on_external_resources = 0.0
    #   Calculate profit expectation
    profit_expectation = sum([lp.profit_npv_per_ha * lp.area_ha for lp in agent.land_parcels])
    #   Calculate profit to date and years in land use for each land parcel
    for land_parcel in agent.land_parcels
        agent.profit_to_date += land_parcel.total_yearly_profit * timedelta
        land_parcel.years_in_current_land_use += timedelta
    end

    # retire on reaching 65, negative profit expectation over whole farm conditional on likelihood to change, or immediately on expectation below 2000/ha average over whole farm
    if (agent.age>=65) || ((profit_expectation < 0) && (agent.likelihood_to_change > rand())) || (profit_expectation / sum([lp.area_ha for lp in agent.land_parcels]) < -2000)
        #   Retire agent and return the new agent
        agent = retire(agent,land_uses,utility_barrier_inputs)
    end
    return agent
end


function step!(agents::Vector{Agent},land_uses::Vector{LandUse},external_resources::Vector{ExternalResource},catchments::Vector{Catchment},timestep,timedelta,years_until_taxes,disable_land_use_change,carbon_price_growth_path,N2OPrice,CH4Price,CO2Price)
    utility_barrier_inputs = DataFrame(XLSX.readtable("data/Model inputs.xlsx","utility_barrier_matrix"))
    agents = advance_time.(agents,Ref(land_uses),Ref(timedelta),Ref(utility_barrier_inputs))
    Learning.learn_by_doing!.(agents,Ref(land_uses))
    Learning.learn_from_neighbours!.(agents,Ref(land_uses))
    Learning.learn_from_resources!(agents,land_uses,external_resources)
    Policies.increase_wq_prices!(agents)
    if carbon_price_growth_path == "exponential"
        Policies.exponential_price_path!(land_uses)
    elseif carbon_price_growth_path == "parabolic"
        Policies.parabolic_price_path!(land_uses,timestep,N2OPrice,CH4Price,CO2Price)
    else
        throw(ArgumentError("Invalid carbon growth path"))
    end
    @time DecisionProblem.manage_land_parcels!.(agents,Ref(land_uses),Ref(disable_land_use_change),Ref(years_until_taxes))
    GenerateOutcomes.agent_and_land_outcomes!.(agents,Ref(land_uses),Ref(years_until_taxes))
    GenerateOutcomes.catchment_outcomes!.(catchments,Ref(land_uses))
    return agents
end

end