module GenerateOutcomes
export calculate_outcomes!
using ..Structs, ..DecisionProblem


function agent_and_land_outcomes!(agent::Agent,land_uses::Vector{LandUse})
    time_horizon = 50
    switching_cost = 0
    for land_parcel in agent.land_parcels
        land_use = filter(x -> x.name == land_parcel.land_use_name, land_uses)[1]
        intensity = land_parcel.intensity
        nl_mitigation = land_parcel.nl_mitigation
        production_capability = land_use.max_production_capability * (( agent.plant_and_animal_capability_by_land_use[land_use.name] ^ 0.5 * agent.people_capability ^ 0.25 * agent.business_capability ^ 0.25) + land_use.production_capability_transform)
        production_capability = clamp(production_capability,1e-3,Inf)
        year = round(Integer,land_parcel.years_in_current_land_use)

        n_mitigation, p_mitigation = DecisionProblem.split_nl_mitigation(nl_mitigation,land_parcel::LandParcel) 
        nitrogen_loss = DecisionProblem.calculate_nutrient_loss(land_use,land_parcel,agent,production_capability,intensity,n_mitigation,"Nitrogen")
        phosphorous_loss = DecisionProblem.calculate_nutrient_loss(land_use,land_parcel,agent,production_capability,intensity,p_mitigation,"Phosphorous")
        
        if land_parcel.uid == "land_parcel_86"
            println("intensity: ",intensity)
            println("sediment: ", land_parcel.sediment_loss)
        end
        sediment_loss = DecisionProblem.calculate_sediment_loss(land_use,land_parcel,intensity)
        if land_parcel.uid == "land_parcel_86"
            println("sediment: ", sediment_loss)
        end

        methane_emissions = DecisionProblem.calculate_methane_emissions(land_use,land_parcel,agent,production_capability,intensity,nl_mitigation)
        nitrous_oxide_emissions = DecisionProblem.calculate_nitrous_oxide_emissions(land_use,land_parcel,agent,intensity,nl_mitigation)

        wq_pollutant_costs, ghg_emission_costs = DecisionProblem.calculate_pollutant_costs(land_use,land_parcel,agent,nitrogen_loss,phosphorous_loss,sediment_loss,methane_emissions,nitrous_oxide_emissions,time_horizon)
            # take only the first element of the returned arrays, which are the costs for the current year
        land_parcel.wq_pollutant_costs = wq_pollutant_costs[1] * land_parcel.area_ha
        if ! (occursin("forestry",land_parcel.land_use_name) || occursin("fallow",land_parcel.land_use_name))

            land_parcel.ghg_emission_costs = ghg_emission_costs[1] * land_parcel.area_ha
            land_parcel.sequestration_payments = 0
            land_parcel.carbon_sequestered = 0
        else
            sequestration_payments = DecisionProblem.calculate_sequestration_payments(land_use,land_parcel,agent,production_capability,intensity,nl_mitigation,year,time_horizon)[1]
            land_parcel.sequestration_payments = sequestration_payments * land_parcel.area_ha
            land_parcel.carbon_sequestered = land_parcel.area_ha * land_parcel.carbon_sequestration_adjustment_by_land_use[land_use.name] * land_use.carbon_sequestration_over_time[year+1]
            land_parcel.ghg_emission_costs = 0
        end
        _, profit_anpv = DecisionProblem.calculate_profit_anpv(land_use,land_parcel,agent,production_capability,intensity,nl_mitigation,year,time_horizon,switching_cost,nitrogen_loss,phosphorous_loss,sediment_loss,methane_emissions,nitrous_oxide_emissions)
        utility = DecisionProblem.calculate_utility(land_use,land_parcel,agent,production_capability,agent.nl_mitigation_capability,intensity,nl_mitigation)
        land_parcel.nitrogen_loss = nitrogen_loss
        land_parcel.phosphorous_loss = phosphorous_loss
        land_parcel.methane_emissions = methane_emissions
        land_parcel.nitrous_oxide_emissions = nitrous_oxide_emissions
        land_parcel.sediment_loss = sediment_loss
         
        land_parcel.profit_npv_per_ha = profit_anpv
        land_parcel.utility = utility

        revenue_at_maturity = DecisionProblem.calculate_revenue(land_use,land_parcel,agent,production_capability,intensity,nl_mitigation)
        costs_at_maturity = DecisionProblem.calculate_variable_costs(land_use,land_parcel,agent,production_capability,intensity,nl_mitigation)
        land_parcel.total_yearly_profit = land_parcel.area_ha * DecisionProblem.calculate_profit_over_time(revenue_at_maturity,costs_at_maturity,land_use.fixed_costs,land_use.years_to_yield,land_use.crop_lifetime)[round(Integer,land_parcel.years_in_current_land_use+1)]
        land_parcel.total_yearly_profit += land_parcel.sequestration_payments
        land_parcel.total_yearly_profit -= land_parcel.wq_pollutant_costs
        land_parcel.total_yearly_profit -= land_parcel.ghg_emission_costs
    end
end


function calculate_out_of_model_sediment_load!(catchment::Catchment,current_period_load_in_model)
    # reduce the out-of-model sediment load by the same proportion as the in-model load relative to the minimum possible load
    if catchment.Load_Sediment > 0
        min_in_model_load = sum([lp.sediment_loss_potenial_by_land_use["native-forestry"] * lp.area_ha for lp in catchment.land_parcels])
        last_period_load_in_model = catchment.Load_Sediment - catchment.out_of_model_Sediment_load
        # avoid floating point precision weirdness
        if (abs(last_period_load_in_model - current_period_load_in_model) < 1e-6) || (abs(last_period_load_in_model - min_in_model_load) < 1e-6)
            return
        end
        proportional_progress = (last_period_load_in_model - current_period_load_in_model) / (last_period_load_in_model - min_in_model_load)
        oom_equivalent_reduction = (catchment.out_of_model_Sediment_load - catchment.min_out_of_model_Sediment_load) * proportional_progress
        catchment.out_of_model_Sediment_load = catchment.out_of_model_Sediment_load - oom_equivalent_reduction
    end
end

function catchment_outcomes!(catchment::Catchment,land_uses::Vector{LandUse})
    total_N_leached = 0.0
    total_P_leached = 0.0
    total_Sediment_leached = 0.0
    prod_area = 0.0
    f_and_f_area = 0.0

    N_leaching_scalar = 0.64 # calculated from mean N leaching in MEIHAWAI output / mean N yield in Snelder et al, to account for attenuation
    P_leaching_scalar = 0.54 # calculated from mean P leaching in MEIHAWAI output / mean P yield in Snelder et al, to account for attenuation
    Sed_leaching_scalar = 1.0 # sediment yields calculated from same underlying model so no need to scale
    baseline_N_leaching_per_ha = sum([lp.N_intensity_coefficients_by_land_use["native-forestry"][1] * lp.area_ha for lp in catchment.land_parcels]) / sum([lp.area_ha for lp in catchment.land_parcels])
    baseline_P_leaching_per_ha = sum([lp.P_intensity_coefficients_by_land_use["native-forestry"][1] * lp.area_ha for lp in catchment.land_parcels]) / sum([lp.area_ha for lp in catchment.land_parcels])
    # sediment base leaching weighted higher than pure forest to account for some bare land making excess contributions
    

    area_by_lu = Dict(lu.name => 0.0 for lu in land_uses)
    for lp in catchment.land_parcels
        total_N_leached += lp.nitrogen_loss * lp.area_ha
        total_P_leached += lp.phosphorous_loss * lp.area_ha
        total_Sediment_leached += lp.sediment_loss * lp.area_ha
        if !(lp.land_use_name in ["fallow","production-forestry","carbon-forestry","native-forestry"])
            prod_area += lp.area_ha
        else
            f_and_f_area += lp.area_ha
        end
        area_by_lu[lp.land_use_name] += lp.area_ha
    end
    total_N_leached += baseline_N_leaching_per_ha * (catchment.watershed_area_ha - (prod_area + f_and_f_area))
    total_P_leached += baseline_P_leaching_per_ha * (catchment.watershed_area_ha - (prod_area + f_and_f_area))

    calculate_out_of_model_sediment_load!(catchment,total_Sediment_leached)
    total_Sediment_leached += catchment.out_of_model_Sediment_load

    catchment.area_in_forestry_or_fallow = f_and_f_area
    catchment.area_in_production = prod_area
    catchment.primary_land_use = sort(collect(area_by_lu),by=x->x[2],rev=true)[1][1]
    catchment.total_coverage = sum([lp.area_ha for lp in catchment.land_parcels]) / catchment.watershed_area_ha
    catchment.Load_N = total_N_leached * N_leaching_scalar
    catchment.Load_P = total_P_leached * P_leaching_scalar
    catchment.Load_Sediment = total_Sediment_leached * Sed_leaching_scalar
    catchment.excess_N_percent = 100 * max(0,(catchment.Load_N - catchment.MAL_N)) / catchment.Load_N
    catchment.excess_P_percent = 100 * max(0,(catchment.Load_P - catchment.MAL_P)) / catchment.Load_P
    catchment.excess_Sediment_percent = 100 * max(0,(catchment.Load_Sediment - catchment.MAL_Sediment)) / catchment.Load_Sediment

    # calculate nutrient loss mitigation subsidies for next timestep
    catchment.tax_revenue_pool = sum([lp.wq_pollutant_costs for lp in catchment.land_parcels])
    for lp in catchment.land_parcels
        if !(lp.land_use_name in ["fallow","production-forestry","carbon-forestry","native-forestry"])
            # subsidy is a per-hectare value, all productive hectares get an equal share of the pool as a subsidy
            lp.nl_mitigation_subsidy = catchment.tax_revenue_pool / catchment.area_in_production
        else
            lp.nl_mitigation_subsidy = 0
        end
    end
end

end