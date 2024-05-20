module DecisionProblem
using ..Structs, ..Learning, Optim, Distributions, Memoize, Statistics
export choose_land_use!, calculate_utility, calculate_environmental_damgage, calculate_profitability

# heterogeneity of response to mitigation 
# heterogeneity of initial outcome/leaching levels
# targeted vs specific = interventions to convince people to do higher mitigation where it's more effective
# intervention - levy - communalise cost of mitigation but target on efficacy   

function setup_optimizer_params(agent::Agent,land_parcel::LandParcel,land_use::LandUse,production_capability)
    
    inbounds = 1e-3

    # analytic solution to 1-d profit optimisation used to set initial condition for 2-d optimiser
    if occursin("forestry",land_use.name) || land_use.name == "fallow"
        # forestry has no variable costs so intensity is irrelevant, set to 1 to ensure max production for production forests
        intensity_ansatz = 1.0
    else
        intensity_ansatz = log(land_parcel.max_revenue_by_land_use[land_use.name] * production_capability / land_use.intensity_variable_cost) / production_capability
    end
    x0 = [intensity_ansatz,0.1]
    
    #   agent.likelihood to change lower and upper bounds
    lower = [land_use.intensity_limits[1],land_use.nl_mitigation_limits[1]]
    upper = [land_use.intensity_limits[2],land_use.nl_mitigation_limits[2]]
    
    #   If agent is non-compliant, reset upper and lower bounds to ignore restrictions
    if agent.is_non_compliant
        upper[1] = 1.0
        lower[2] = 0.0
    end
    
    #   Check if any x0 parameters are beyond lower or upper bounds, if so, adjust to be within bounds
    for (i,param) in enumerate(x0)
        if param < lower[i]
            x0[i] = lower[i] + inbounds
        end
        if param > upper[i]
            x0[i] = upper[i] - inbounds
        end
    end
    return x0, lower, upper
end

function calculate_revenue(land_use::LandUse,land_parcel::LandParcel,agent::Agent,production_capability,intensity,nl_mitigation)
    #   Calculates max_revenue (attribute of land parcel and specific land_use), adjusts considers capacility and intensity
    max_revenue = land_parcel.max_revenue_by_land_use[land_use.name] * land_use.revenue_adjustment
    return max_revenue * (1-exp(-production_capability * intensity))
end


function calculate_variable_costs(land_use::LandUse,land_parcel::LandParcel,agent::Agent,production_capability,intensity,nl_mitigation)
    #   Adjusts max variable costs for land use and considers intensity and non-linear mitigation
    return intensity * land_use.intensity_variable_cost + max(0, nl_mitigation * land_use.nl_mitigation_variable_cost - land_parcel.nl_mitigation_subsidy)
end

function calculate_nutrient_loss(land_use::LandUse,land_parcel::LandParcel,agent::Agent,nl_mitigation_capability,intensity,nl_mitigation,nutrient)
    
    #   max_nl = calculated for nutrient type (Nitrogen, Phosphorous) - attribute of speicifc land_use for given land_parcel
    
    a,b,c = nutrient == "Nitrogen" ? land_parcel.N_intensity_coefficients_by_land_use[land_use.name] : land_parcel.P_intensity_coefficients_by_land_use[land_use.name]
    mitigation_potential = nutrient == "Nitrogen" ? land_parcel.n_mitigation_potential_by_land_use[land_use.name] : land_parcel.p_mitigation_potential_by_land_use[land_use.name]
    # restrict nutrient loss to be within the bounds of the land use
    lower_limit = a + b * land_parcel.min_intensity_by_land_use[land_use.name] + c * land_parcel.min_intensity_by_land_use[land_use.name]^2
    upper_limit = a + b * land_parcel.max_intensity_by_land_use[land_use.name] + c * land_parcel.max_intensity_by_land_use[land_use.name]^2
    #   non-linear nutrient loss - constant and linear coefficients
    #   Nutrient loss = max_nl, nl_mitigation, intensity, nl_coefficients
    return (1 - mitigation_potential * (1-exp(-nl_mitigation_capability*nl_mitigation))) * clamp(a + b*intensity + c * intensity^2,lower_limit,upper_limit)
end


function calculate_methane_emissions(land_use::LandUse,land_parcel::LandParcel,agent::Agent,production_capability,intensity,nl_mitigation)
    #   Calculates methane emissions considering underlying rate, production capability, and intensity
    max_revenue = land_parcel.max_revenue_by_land_use[land_use.name] * land_use.revenue_adjustment
    return land_use.methane_emission_rate * max_revenue * (1-exp(-production_capability * intensity))
end

function calculate_nitrous_oxide_emissions(land_use::LandUse,land_parcel::LandParcel,agent::Agent,intensity,nl_mitigation)
    a,b,c = land_parcel.N2O_intensity_coefficients_by_land_use[land_use.name]
    #   Calculates nitrous oxide emissions considering underlying max rate, intensity, and non-linear coefficients
    lower_limit = a + b * land_parcel.min_intensity_by_land_use[land_use.name] + c * land_parcel.min_intensity_by_land_use[land_use.name]^2
    upper_limit = a + b * land_parcel.max_intensity_by_land_use[land_use.name] + c * land_parcel.max_intensity_by_land_use[land_use.name]^2
    return clamp(a + b*intensity + c * intensity^2,lower_limit,upper_limit)
end

function calculate_sediment_loss(land_use::LandUse,land_parcel::LandParcel,intensity)
    #   Calculates sediment loss considering underlying rate, intensity, and land use
    a,b = land_parcel.Sediment_intensity_coefficients_by_land_use[land_use.name]
    lower_limit = a + b * land_parcel.min_intensity_by_land_use[land_use.name]
    upper_limit = a + b * land_parcel.max_intensity_by_land_use[land_use.name]
    return clamp(a + b*intensity,lower_limit,upper_limit) * land_parcel.sediment_loss_potenial_by_land_use[land_use.name]
end

function calculate_profit_over_time(revenue_at_maturity,costs_at_maturity,fixed_costs,years_to_yield,crop_lifetime)
    #   If Lifetime = 0:
    if crop_lifetime == 0
        #   profit is revenue_matruity - var_cost_maturity - fixed cost
        profit = revenue_at_maturity - costs_at_maturity - fixed_costs
        #   Repeat profit 200 times (200 years)
        profit = repeat([profit], 200)
        return profit
    end
    
    #   If years_to_yield = crop_lifetime (special case for production forestry with a yield every 28 years)
    if crop_lifetime == years_to_yield
        #   Revenue & variable costs during lifetime = 0. On end concatenate revenue and variable costs for maturity year 
        revenue = vcat(zeros(crop_lifetime), [revenue_at_maturity])
        var_costs = vcat(zeros(crop_lifetime), [costs_at_maturity])

        #   Fixed costs incurred in all years of crop lifetime
        fixed_costs = fixed_costs * ones(crop_lifetime+1)

        profit = revenue - var_costs - fixed_costs

        #   Repeats profit cycle over crop lifetime for number of cycles in 200 years
        profit = repeat(profit, (200 ÷ (crop_lifetime)))
        
        return profit
    end

    #   If crop lifetime /= years to yield:
    #   Array 1: length = years_yield, values = 0 -> maturity revenue/ (cost) (Linearly increase)
    #   Array 2: length = lifetime - years_yield, values = maturity revenue/ (cost)
    #   Concatenate arrays
    revenue = vcat(LinRange(0,revenue_at_maturity,years_to_yield+1), revenue_at_maturity * ones(crop_lifetime-years_to_yield))
    var_costs = vcat(LinRange(0,costs_at_maturity,years_to_yield+1), costs_at_maturity * ones(crop_lifetime-years_to_yield))
    
    #   Fixed costs are incurred in all years
    fixed_costs = fixed_costs * ones(crop_lifetime+1)
   
    profit = revenue - var_costs - fixed_costs

    #   Repeats profit cycle over crop lifetime for number of cycles in 100 years
    profit = repeat(profit, (200 ÷ (crop_lifetime)))
    return profit
end

@memoize function get_price_growth_array(growth_factor::Float64,time_horizon)
    #   Calculates price growth array for given growth rate and time horizon
    price_growth_array = [growth_factor^(i-1) for i in 1:time_horizon]
    return price_growth_array
end

function calculate_pollutant_costs(land_use::LandUse,land_parcel::LandParcel,agent::Agent,nitrogen_loss,phosphorous_loss,sediment_loss,methane_emissions,nitrous_oxide_emissions,time_horizon,years_until_taxes)

    wq_growth_rate = 1.03
    wq_price_growth_array = get_price_growth_array(wq_growth_rate,time_horizon)
    ghg_price_growth_array = get_price_growth_array(1 + agent.carbon_price_growth_rate_forecast,time_horizon)
   
    ghg_tax = land_use.methane_price * methane_emissions + land_use.nitrous_oxide_price * nitrous_oxide_emissions

    nitrogen_tax = nitrogen_loss * land_parcel.nitrogen_price
    phosphorous_tax = phosphorous_loss * land_parcel.phosphorous_price
    sediment_tax = sediment_loss * land_parcel.sediment_price

    # no WQ tax if already at best possible state for that contaminant, i.e., forest or fallow for N and P, forestry only for sediment.
    if occursin("forestry",land_use.name) || occursin("fallow",land_use.name)
        nitrogen_tax = 0
        phosphorous_tax = 0
        if occursin("forestry",land_use.name)
            sediment_tax = 0
        end
    end
    wq_tax = nitrogen_tax + phosphorous_tax + sediment_tax

    wq_costs_over_time = wq_tax * wq_price_growth_array
    ghg_costs_over_time = ghg_tax * ghg_price_growth_array

    # zero out wq and ghg taxes until they take effect
    wq_costs_over_time[1:years_until_taxes] .= 0
    ghg_costs_over_time[1:years_until_taxes] .= 0
    # #   Taxes adjusted for growth rate pushed to tax over time

    #   Discount pollutant costs by 1/3 for non-compliant agents
    if agent.is_non_compliant
        wq_costs_over_time .*= 2/3
        ghg_costs_over_time .*= 2/3
    end
    #   Discounted costs over time = costs over time / (1+r)^index
    financial_discount_array = get_price_growth_array(1/(1 + agent.financial_discount_rate),time_horizon)
    wq_discounted = wq_costs_over_time .* financial_discount_array
    ghg_discounted = ghg_costs_over_time .* financial_discount_array
    return wq_discounted, ghg_discounted
end

function calculate_sequestration_payments(land_use::LandUse,land_parcel::LandParcel,agent::Agent,production_capability,intensity,nl_mitigation,years_in_current_land_use,time_horizon)
    #   Initialise empty array for sequestration payments over time
    
    carbon_price_growth_array = get_price_growth_array(1 + agent.carbon_price_growth_rate_forecast,time_horizon)
    
    sequestration_over_time = land_use.carbon_sequestration_over_time[years_in_current_land_use+1:years_in_current_land_use+time_horizon]
    adjusted = sequestration_over_time * land_parcel.carbon_sequestration_adjustment_by_land_use[land_use.name] * 1000 * land_use.carbon_price .* carbon_price_growth_array
    #   For each year:
    #   Discount the sequestration payments to NPV
    discount_array = get_price_growth_array(1/(1 + agent.financial_discount_rate),time_horizon)
    discounted = adjusted .* discount_array
    return discounted
end

function calculate_profit_anpv(land_use,land_parcel,agent,production_capability,intensity,nl_mitigation,years_in_current_land_use,time_horizon,switching_cost,nitrogen_loss,phosphorous_loss,sediment_loss,methane_emissions,nitrous_oxide_emissions,years_until_taxes)
    #   Calculates profit. Considers land use, land parcel, agent, production capability, intensity, non-linear mitigation, years in current land use, time horizon, switching cost, nitrogen loss, phosphorous loss, methane emissions, and nitrous oxide emissions.
    #   Returns revenue at maturity and variable costs at maturity
    revenue_at_maturity = calculate_revenue(land_use,land_parcel,agent,production_capability,intensity,nl_mitigation)
    costs_at_maturity = calculate_variable_costs(land_use,land_parcel,agent,production_capability,intensity,nl_mitigation)
    
    #   Calculate profitability over time for 100 years for given crop situation
    profitability_over_time = calculate_profit_over_time(revenue_at_maturity,costs_at_maturity,land_use.fixed_costs,land_use.years_to_yield,land_use.crop_lifetime)
   
    #   Indexed tuples of profit in years of time horizion, PV profit = FV/(1+r)^index, cumulative PV of profits over time horizion calculated
    #   NPV = cumulative PV of profits - switching cost (opportunity costs of switching land use)
    discount_rate_array = get_price_growth_array(1/(1 + agent.financial_discount_rate),time_horizon)
    profit_npv = sum(profitability_over_time[years_in_current_land_use+1:years_in_current_land_use+time_horizon] .* discount_rate_array) - switching_cost
    
   #   If land use is forestry, account for carbon sequestration benefit

    #   Otherwise, account for pollutant costs
    wq_pollutant_costs_over_time, ghg_emission_costs_over_time = calculate_pollutant_costs(land_use,land_parcel,agent,nitrogen_loss,phosphorous_loss,sediment_loss,methane_emissions,nitrous_oxide_emissions,time_horizon,years_until_taxes)

    profit_npv -= sum(wq_pollutant_costs_over_time)
    if occursin("forestry",land_use.name)
        #   Calculate PV of sequestration payments and add to profit_npv
        sequestration_payments_over_time = calculate_sequestration_payments(land_use,land_parcel,agent,production_capability,intensity,nl_mitigation,years_in_current_land_use,time_horizon)
        profit_npv += sum(sequestration_payments_over_time)
    else
        profit_npv -= sum(ghg_emission_costs_over_time)
    end
    #    Convert profit NPV to annualised profit using an annuity scalar
    annuity_scalar = 1/sum(discount_rate_array)
    annualised_profit = annuity_scalar * profit_npv
    return annuity_scalar, annualised_profit
end


function calculate_environmental_npv(weighting,emissions,annualised_profit,time_horizon,discount_rate)
    #   Calculates environmental npv. Considers weighting, emissions, annualised profit, time horizon, and discount rate.
    profit_scalar = clamp(annualised_profit,0,Inf)
    discount_vector = get_price_growth_array(1/(1+discount_rate),time_horizon)
    return sum(weighting * profit_scalar * emissions * discount_vector)
end

function split_nl_mitigation(nl_mitigation,land_parcel::LandParcel)
    # if no prices set, apply same level of mitigation to each nutrient. 
    # otherwise, mitigate each nutrient in proportion to its past year's contribution to costs, capped at mitigation = 1.0
    total_cost = land_parcel.nitrogen_price * land_parcel.nitrogen_loss + land_parcel.phosphorous_price * land_parcel.phosphorous_loss
    if total_cost == 0
        return nl_mitigation, nl_mitigation
    end
    n_proportion  = land_parcel.nitrogen_price * land_parcel.nitrogen_loss / total_cost
    p_proportion  = land_parcel.phosphorous_price * land_parcel.phosphorous_loss / total_cost
    n_mitigation = 2*n_proportion * nl_mitigation
    p_mitigation = 2*p_proportion * nl_mitigation
    # if either mitigation quantity exceeds 1.0, redistribute the excess to the other nutrient
    if n_mitigation > 1.0
        p_mitigation += (n_mitigation - 1.0) 
        n_mitigation = 1.0
    end
    if p_mitigation > 1.0
        n_mitigation += (p_mitigation - 1.0) 
        p_mitigation = 1.0
    end
    return n_mitigation, p_mitigation
end

function calculate_utility(land_use::LandUse,land_parcel::LandParcel,agent::Agent,production_capability,nl_mitigation_capability,intensity,nl_mitigation,years_until_taxes)
    #   Calculates utility of given land parcel for given agent
    #   Considers land use, production capability, nl_mitigation capability, and intensity

    time_horizon = 50
    #   Gets switching cost of land use for specifc land_parcel - stored in dict attribute of land_use object
    switching_cost = land_use.cost_to_switch_from[land_parcel.land_use_name]

    #   years land_parcel in current use calculated:
    #   If in best use => number of years parcel in current use, otherwise => 0
    years_in_current_land_use = land_parcel.land_use_name == land_use.name ? round(Integer,land_parcel.years_in_current_land_use) : 0
    
    #   Nutrient loss of agents land parcels with given land uses returned
    nitrogen_mitigation, phosphorous_mitigation = split_nl_mitigation(nl_mitigation,land_parcel)
    nitrogen_loss = calculate_nutrient_loss(land_use,land_parcel,agent,nl_mitigation_capability,intensity,nitrogen_mitigation,"Nitrogen")
    phosphorous_loss = calculate_nutrient_loss(land_use,land_parcel,agent,nl_mitigation_capability,intensity,phosphorous_mitigation,"Phosphorous")
    sediment_loss = calculate_sediment_loss(land_use,land_parcel,intensity)
    
    #   Methane and nitrous oxide emissions of agents land parcels with given land uses returned
    methane_emissions = calculate_methane_emissions(land_use,land_parcel,agent,production_capability,intensity,nl_mitigation)
    nitrous_oxide_emissions = calculate_nitrous_oxide_emissions(land_use,land_parcel,agent,intensity,nl_mitigation)


    #   Returns annualised profit and environmental npv for given land use, land parcel, agent, production capability, nl_mitigation capability, intensity, non-linear mitigation, years in current land use, time horizon, switching cost, nitrogen loss, phosphorous loss, methane emissions, and nitrous oxide emissions.
    annuity_scalar, annualised_profit = calculate_profit_anpv(land_use,land_parcel,agent,production_capability,intensity,nl_mitigation,years_in_current_land_use,time_horizon,switching_cost,nitrogen_loss,phosphorous_loss,sediment_loss,methane_emissions,nitrous_oxide_emissions,years_until_taxes)
    
    #   Calculates npv of environmental losses. Considers weighting, emissions, annualised profit, time horizon, and discount rate.
    nitrogen_loss_npv, phosphorous_loss_npv, methane_emissions_npv, nitrous_oxide_emissions_npv, sediment_loss_npv = calculate_environmental_npv.([agent.nitrogen_weighting,agent.phosphorous_weighting,agent.ghg_weighting,agent.ghg_weighting,agent.sediment_weighting],
                                                                                                                                                                    [nitrogen_loss,phosphorous_loss,methane_emissions,nitrous_oxide_emissions,sediment_loss],
                                                                                                                                                                  Ref(annualised_profit),Ref(time_horizon),Ref(agent.env_discount_rate))
    #   Utility = annualised profit - annuity scalar * (sum of environmental npvs), environmental NPVS * scalar to get annual environmental cost
    utility = annualised_profit - annuity_scalar * (nitrogen_loss_npv + phosphorous_loss_npv + methane_emissions_npv + nitrous_oxide_emissions_npv + sediment_loss_npv)
    return utility
end

function calculate_production_capability(agent,land_use,forecast_learning_adjustment;time_horizon=50)

    #   Stores the pre-adjustment values of the agent's capabilities in the preadjustment vector.
    preadjustment = [agent.plant_and_animal_capability_by_land_use[land_use.name],agent.people_capability,agent.business_capability,agent.nl_mitigation_capability]
    
    #   If forecast_learning = true:
    if forecast_learning_adjustment

        #   inv_dr = value used to calculate PV of FCF and NPV
        inv_dr = 1/(1+agent.financial_discount_rate)
        
        #   calculates the midpoint of the time horizon
        midpoint = log(inv_dr,0.5*inv_dr^(time_horizon+1)+0.5)
        midpoint = round(Integer,midpoint)

        #   iterates over the midpoint
            for i in 1:midpoint
                #   updates the agent's capabilities based on the agents activity at each point in time
                Learning.learn_by_doing!(agent,[land_use])
            end
    end

    #   Weighted sum of the agent's capabilities, with different weights for people and business capabilities.
    production_capability = agent.plant_and_animal_capability_by_land_use[land_use.name] ^ 0.5 * agent.people_capability ^ 0.25 * agent.business_capability ^ 0.25
    production_capability = clamp(production_capability+land_use.production_capability_transform,1e-3,Inf)
    #   Capability adjusted on the agent's deviation parameter and the land use's maximum production capability.
    production_capability *= agent.deviation_parameter * land_use.max_production_capability

    #   Agent's capabilities are reset to their pre-adjustment values, and the function returns the production_capability
    agent.plant_and_animal_capability_by_land_use[land_use.name] = preadjustment[1]
    agent.people_capability = preadjustment[2]
    agent.business_capability = preadjustment[3]
    agent.nl_mitigation_capability = preadjustment[4]

    
    return production_capability
end

function value_land_use(agent::Agent,land_parcel::LandParcel,land_use::LandUse,years_until_taxes;forecast_learning_adjustment=false,tolerance=5e-3)
#   optimization technique: optimal intensity and nl_mitigation that maximises utility

    #   Calculate production capability - weighted sum of people and business capability and considers activity learning
    production_capability = calculate_production_capability(agent,land_use,forecast_learning_adjustment)
    
    #   Calculate non-linear mitigation capability and adjust for deviation parameter
    nl_mitigation_capability = agent.nl_mitigation_capability * land_use.max_nl_mitigation_capability * agent.deviation_parameter
    
    function optimise_value(params::Vector{Float64})
    #   Nested function - optimise_value: Calculates utility for intensity and nl_mitigation    
        intensity = params[1]
        nl_mitigation = params[2]

        #   If intensity or nl_mitigation invalid, return infinity
        if intensity < 0 || intensity > 1 || nl_mitigation < 0 || nl_mitigation > 1
            return Inf
        end
        #   Otherwise, calculate utility with function calculate_utility
        utility = calculate_utility(land_use,land_parcel,agent,production_capability,nl_mitigation_capability,intensity,nl_mitigation,years_until_taxes)
        return -1 * utility
    end

    #   Set up optimizer parameters - x0, lower, upper
    x0, lower, upper = setup_optimizer_params(agent,land_parcel,land_use,production_capability)
    local intensity
    local nl_mitigation
    try
        #   Optimise functiones takes function to be minimised, lower and upper bounds, initial conditions, and optimization technique
        res = optimize(optimise_value,
                        lower,
                        upper,
                        x0,
                        Fminbox(GradientDescent()),
                        Optim.Options(f_calls_limit=250,
                                    outer_iterations=10,
                                    f_tol=tolerance,
                                    g_tol=tolerance,
                                    x_tol=tolerance,))
        intensity, nl_mitigation = res.minimizer
    catch
        #   If optimisation fails, return initial conditions
        println("OPTIMIZER FAILED FOR ", land_parcel.land_use_name, " to ", land_use.name, " with initial conditions ", x0, " and production capability ", production_capability)
        println(land_parcel.years_in_current_land_use)
        println(land_parcel.land_use_name)
        intensity = x0[1]
        nl_mitigation = x0[2]
    end
    if intensity > land_use.intensity_limits[2]
        intensity = land_use.intensity_limits[2] + (intensity - land_use.intensity_limits[2]) * 1/3
    end
    if nl_mitigation < land_use.nl_mitigation_limits[1]
        nl_mitigation = land_use.nl_mitigation_limits[1] - (land_use.nl_mitigation_limits[1] - nl_mitigation) * 1/3
    end

    #   Recall calculate_utility function with intensity and nl_mitigation updated 
    utility = calculate_utility(land_use,land_parcel,agent,production_capability,nl_mitigation_capability,intensity,nl_mitigation,years_until_taxes)

    return utility, intensity, nl_mitigation
end


function restrict_change(land_parcel,proportion_of_change,intensity,nl_mitigation)
    #   Set intensity and nl_mitigation to land_parcel's attribute values + adjustment for potential changes
    intensity = land_parcel.intensity + proportion_of_change * (intensity - land_parcel.intensity)
    nl_mitigation = land_parcel.nl_mitigation + proportion_of_change * (nl_mitigation - land_parcel.nl_mitigation)
    return intensity, nl_mitigation
end

function manage_existing_land_use!(agent::Agent,land_parcel::LandParcel,land_uses::Vector{LandUse},proportion_of_change::Float64,years_until_taxes)
    #   Manages existing land use of land parcel. Takes agent, land parcel, land uses, and initialise as arguments. Calculates utility, intensity, and nl_mitigation for existing use
    util, intensity, nl_mitigation, = value_land_use(agent,land_parcel,filter(x->x.name==land_parcel.land_use_name,land_uses)[1],years_until_taxes,tolerance=1e-5)
        #   Otherwise, maintain as attribute values + change potential adjustment
    intensity,nl_mitigation = restrict_change(land_parcel,proportion_of_change,intensity,nl_mitigation)
    #   Set attributes of land parcel to calculated values
    land_parcel.intensity = intensity
    land_parcel.nl_mitigation = nl_mitigation
end

function evaluate_land_uses!(agent::Agent,land_parcel::LandParcel,land_uses::Vector{LandUse},proportion_of_change::Float64,disallow_forestry,years_until_taxes)
    #   Evaluating different land uses for a given land parcel and deciding on the most beneficial one. 
    
    #   Checks if agent is initalising. If true, years_in_lu_adjustment = random number 1 -> 18, false, years_in_lu_adjustment = 0
    years_in_lu_adjustment = disallow_forestry ? rand(1:round(Integer,agent.age-18)) : 0
    #   Initalising = true, filter out unallowed land uses and forestry for land_parcecl 
    #   Initalising = false, filter out unallowed land uses for land_parcel
    allowed_land_uses = disallow_forestry ? filter(lu -> (lu.name in land_parcel.allowed_land_uses) & (! occursin("forestry",lu.name)) & (lu.name != "fallow") ,land_uses) : filter(lu->lu.name in land_parcel.allowed_land_uses,land_uses)
    
    # handles the edge case where allowed land uses ends up empty (e.g., pixels on water's edge which have no economic indicator values)
    if isempty(allowed_land_uses)
        allowed_land_uses = filter(lu->lu.name == "fallow",land_uses)
    end
    #   Calculate land_use_values, intensities, and nl_mitigations for each allowed land use -  value_land_use function - returns optimum utility, intensity, and nl_mitigation
    prospective_land_use_values, intensities, nl_mitigations = zip(value_land_use.(Ref(agent),Ref(land_parcel), allowed_land_uses, Ref(years_until_taxes), forecast_learning_adjustment=true)...)
    
    # Finds the index of the current land use in the list of allowed land uses for land_parcel
    current_land_use_idx = findfirst(x->x==land_parcel.land_use_name,[lu.name for lu in allowed_land_uses])

    current_land_use_category = filter(x->x.name==land_parcel.land_use_name,land_uses)[1].category
    
    utility_adjustments = [agent.utility_barrier_matrix[current_land_use_category][filter(x->x.name==lu.name,land_uses)[1].category] for lu in allowed_land_uses]
   
    if current_land_use_idx === nothing
        comparative_utilities = [prospective_land_use_values[i] for i in eachindex(prospective_land_use_values)]
    else
        
        #   Otherwise, calculate the difference between utility of prospective_land_use_values and utility of the land_use at current_land_use_idx, scaled by the agent's utility barrier matrix
        comparative_utilities = [(prospective_land_use_values[i] / utility_adjustments[i]) - prospective_land_use_values[current_land_use_idx] for i in eachindex(prospective_land_use_values)]
    end
    if (maximum(comparative_utilities) > 0) | (current_land_use_idx === nothing)
        
        #  Get index of maximum utility, get name of land_use at that index. and get the land_use for land parcel to it
        
        lu_index = argmax(comparative_utilities)
        new_lu = allowed_land_uses[lu_index]
        # no switching costs in first timestep
        if ! agent.undertaking_initial_evaluation
            land_parcel.switching_costs_incurred = new_lu.cost_to_switch_from[land_parcel.land_use_name] * land_parcel.area_ha
        end
        land_parcel.land_use_name = new_lu.name
        
        
        #   If new land_use = forestry, restrict allowed land uses to prevent harvesting/switching out
        if occursin("forestry",new_lu.name)
            land_parcel.allowed_land_uses = [new_lu.name]
        end
        #   Recall manage land_use with updated land_parcel, land_uses, and initialise properties
        land_parcel.years_in_current_land_use = years_in_lu_adjustment
        manage_existing_land_use!(agent,land_parcel,land_uses,1.0,years_until_taxes)

    else
        #   Otherwise, manage existing land use of non-switched land_use
        manage_existing_land_use!(agent,land_parcel,land_uses,proportion_of_change,years_until_taxes)
    end
end


function get_evaluation_probability_vector(agent::Agent,land_uses::Vector{LandUse},years_until_taxes,disable_land_use_change)
    if agent.undertaking_initial_evaluation
        return ones(length(agent.land_parcels))
    elseif disable_land_use_change
        return zeros(length(agent.land_parcels))
    end
    # recalculate utilities using last period's intensity and mitigation 
    current_utilities = []
    for lp in agent.land_parcels
        current_lu = filter(x->x.name==lp.land_use_name,land_uses)[1]
        production_capability = calculate_production_capability(agent,current_lu,false)
        nl_mitigation_capability = agent.nl_mitigation_capability * current_lu.max_nl_mitigation_capability * agent.deviation_parameter
        current_utility = calculate_utility(current_lu,lp,agent,production_capability,nl_mitigation_capability,lp.intensity,lp.nl_mitigation,years_until_taxes)
        push!(current_utilities,current_utility)
    end
    max_util = max(100, 1.15 * maximum(current_utilities)) # 15% buffer to utility ceiling to allow for evaluation of best parcel
    relative_utils = max_util .- current_utilities
    # cap relative likelihood of evaluation at 5x the mean
    prob_vec = relative_utils ./ sum(relative_utils) .* agent.likelihood_to_change .* length(agent.land_parcels)
    # zero out evaluation probability if only one land use is allowed
    mask = [length(lp.allowed_land_uses) == 1 ? 0 : 1 for lp in agent.land_parcels]
    prob_vec = clamp.(prob_vec .* mask,Ref(0),Ref(1))
    return prob_vec
end


function manage_land_parcels!(agent::Agent,land_uses::Vector{LandUse},disable_land_use_change,years_until_taxes;disallow_forestry=false)
    #   Takes an Agent object and a vector of LandUse objects as arguments, to evaluate land_use's intensity & nl_mitigation for each land parcel managed by the agent.
    #   Should agent evaluate land use (make a land use change decision), or manage existing use:
    #   Generate random number, if likelihood to change > number or initialise = true, => evaluate land use
    evaluation_probability_vector = get_evaluation_probability_vector(agent,land_uses,years_until_taxes,disable_land_use_change)
    for (land_parcel,probability) in zip(agent.land_parcels,evaluation_probability_vector)
        land_parcel.switching_costs_incurred = 0
        chooses_to_evaluate = probability > rand()
        if chooses_to_evaluate
        #   Iterate over land_parcels managed by the agent, and evaluate land use of each
            #   Evaluate land use of land parcel in terms of it's utility and nl_mitigation and whether it is the best opportunity
            evaluate_land_uses!(agent,land_parcel,land_uses,probability,disallow_forestry,years_until_taxes)
    #   Otherwise => manage land use and return intensity and nl_mitigation of land parcel
        else
            manage_existing_land_use!(agent,land_parcel,land_uses,probability,years_until_taxes)
        end
    end
    #  Set agent's initial evaluation to false since they've just completed it
    agent.undertaking_initial_evaluation = false
end

end