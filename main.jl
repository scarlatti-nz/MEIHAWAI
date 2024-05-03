using Revise
includet("src/Structs.jl")
includet("src/Learning.jl")
includet("src/SocialNetwork.jl")
includet("src/Population.jl")
includet("src/DecisionProblem.jl")
includet("src/Policies.jl")
includet("src/GenerateOutcomes.jl")
includet("src/ModelStep.jl")
includet("src/RecordData.jl")



using Main.Structs, Main.Population, Main.Policies, Main.ModelStep, Main.DecisionProblem, Main.RecordData, Main.SocialNetwork, Main.Learning, Distributions, Random, DataFrames, CSV, Statistics, Plots, Dates, ArgParse


Random.seed!(1234)

s = ArgParseSettings()
@add_arg_table s begin
    "--enable-bolus" 
        help = "Enable bolus land uses" 
        arg_type=Bool
        default = false
    "--disable-land-use-change"
        help = "Disable land use change"
        arg_type=Bool
        default = false
    "--extension-scale"
        help = "Scale factor for extension spending (0.0 = no govt spending)"
        arg_type=Float64
        default = 0.0
    "--NPrice"
        help = "Price per kg for nitrogen excesses"
        arg_type=Float64
        default = 0.0
    "--PPrice"
        help = "Price per kg for phosphorous excesses"
        arg_type=Float64
        default = 0.0
    "--SedimentPrice"
        help = "Price per ton for sediment excesses"
        arg_type=Float64
        default = 0.0
    "--N2OPrice"
        help = "Price per ton CO2e for N2O emissions"
        arg_type=Float64
        default = 0.0
    "--CH4Price"
        help = "Price per ton CO2e for CH4 emissions"
        arg_type=Float64
        default = 0.0
    "--CO2Price"
        help = "Price per ton CO2 sequestration"
        arg_type=Float64
        default = 70.0
    "--sediment-model"
        help = "Sediment model to use. default assumes land cover factor 0.466 for trees, alt assumes 0.1"
        arg_type=String
        default = "default"
    "--carbon-growth-path"
        help = "Sets the shape of the carbon growth path. Options are 'parabolic' or 'exponential'"
        arg_type=String
        default = "parabolic"
    "--geodata"
        help = "Data directory"
        arg_type=String
        default = "data/"
    "--output-dir"
        help = "Output directory"
        arg_type=String
        default = "output/"
    "--log-file"
        help = "Log file"
        arg_type=String
        default = "/"
    "--max-parcels"
        help = "Maximum number of parcels to run"
        arg_type=Int64
        default = 5000
    "--sensitivity-analysis"
        help = "Sensitivity analysis scenario to run"
        arg_type=Int64
        default = 0
    "--initialise-only"
        help = "Abort model after initialisation timestep"
        arg_type=Bool
        default=false
end
parsed_args = parse_args(ARGS, s)
geodata = parsed_args["geodata"]
output_folder = parsed_args["output-dir"]

function main(;timedelta = 1)
    ##
    agents, land_uses, land_info, catchments = Population.generate_initial_population(geodata,parsed_args["max-parcels"],parsed_args["sediment-model"])
    ##
    SocialNetwork.assign_neighbours!(agents)
    ##
    ##
    DecisionProblem.manage_land_parcels!.(agents,Ref(land_uses),parsed_args["disable-land-use-change"],disallow_forestry=true)

    ModelStep.add_initial_capability!.(agents,Ref(land_uses))

    GenerateOutcomes.agent_and_land_outcomes!.(agents,Ref(land_uses))

    GenerateOutcomes.catchment_outcomes!.(catchments,Ref(land_uses))

    # adjust down likelihood to change for sheep and beef farmers based on observed patterns of land use change.
    # Population.adjust_l2c!(agents)


    Policies.apply_incentives!(land_uses,agents,catchments,parsed_args["NPrice"],parsed_args["PPrice"],parsed_args["SedimentPrice"],parsed_args["N2OPrice"],parsed_args["CH4Price"],parsed_args["CO2Price"])
    
    baseline_external_resources = Policies.generate_baseline_external_resources(agents)
    additional_extenral_resources = Policies.generate_additional_external_resources(agents,parsed_args["extension-scale"])
    external_resources = [baseline_external_resources;additional_extenral_resources]
    #   Specifies how regulations applied depend on rest of model - deprecated

    # Policies.apply_regulations!(land_uses,agents,parsed_args["regulation-scale"])
    
    #   Accounts for bolus land_uses if bolus is enabled
    Policies.enable_bolus!(land_uses,agents,parsed_args["enable-bolus"])

    if parsed_args["sensitivity-analysis"] > 0
        Population.change_params_for_SA!(parsed_args["sensitivity-analysis"],agents,land_uses)
    end

    RecordData.export_agents(agents,land_uses,land_info,0,output_folder)

    RecordData.export_catchments(catchments,0,output_folder)

    RecordData.export_land_use(agents,land_info,output_folder)
    ##
    
    if parsed_args["initialise-only"]
        open(output_folder*"/complete.txt","w") do log
            write(log,"finished")
        end
        return "aborting after initialisation"
    end
    
    for timestep in 1:25

        #   Advance time and perform model calculations at each step and record the data for each agent
        agents = ModelStep.step!(agents,land_uses,external_resources,catchments,timestep,timedelta,parsed_args["disable-land-use-change"],parsed_args["carbon-growth-path"],parsed_args["N2OPrice"],parsed_args["CH4Price"],parsed_args["CO2Price"])
        RecordData.export_agents(agents,land_uses,land_info,timestep,output_folder)
        RecordData.export_catchments(catchments,timestep,output_folder)
        RecordData.export_extension_spend!(external_resources,timestep,output_folder)
    end
end

open(parsed_args["log-file"],"w") do io
    redirect_stderr(io) do
        try
            main()
        finally
            println("flushing")
            flush(io)
            open(output_folder*"/complete.txt","w") do log
                write(log,"finished")
            end
        end
    end
end