module SocialNetwork
using ..Structs, Random, Statistics, StatsBase, LinearAlgebra, Plots, JSON, MetaGraphs, GraphPlot, GeoDataFrames, ArchGDAL, Distributions

# USE soft random geometric graph with exponential function beta * exp(-rij/r0) with normalisation at the national level such that mean number of connections is 10. This will give sparser and denser segments lower and higher mean degrees respectively.
# use r0 to scale network density (higher r0 = lower density) and then normalise beta such that average connections = 10
# the 1/dij thing fails cause you get alphas close to 1 blowing up the average alpha, which isn't a problem in the exponential version. 

function run(agents::Vector{Agent}; n=2.0, r0=1500)

    #   hcat, horizontally connects arrays; array1 = agent_x, array2 = agent_y
    positions = hcat([a.x_coordinate for a in agents], [a.y_coordinate for a in agents])
    
    #   Number of agents
    npeople = length(positions[:,1])

    #   dict: key = int 1-> npeople (agent number), value = initalised as empty (index of adjacent agents). Stores adjacency information.
    adjacency = Dict(
        i => Int64[] for i in 1:npeople
    )
    
    #   For each agent, calculate the distance to all other agents.
    for j in 1:npeople
        #    for each agent, calculate the distance to all other agents (Euclidean formula)
        distances = sqrt.(
                (positions[:, 1] .- positions[j, 1]) .^ 2
                + (positions[:, 2] .- positions[j, 2]) .^ 2
            )
        #   Probability of each agent being a neighbor, r0 = radius of influence
        prob = n * exp.(-(distances ./ r0))

        #   Select agents to be neighbors based on probability threshold
        selected_neighbours = [n for n in j+1:npeople if prob[n] > rand()]

        #   Add neighbour agents to adjacency list and their selected neighbor
        adjacency[j] = vcat(adjacency[j], selected_neighbours)
        for o in selected_neighbours
            adjacency[o] = vcat(adjacency[o], [j])
        end
    end

    #   creates a list of the lengths of adjacency lists (number of neighbours)
    degrees = [length(adjacency[j]) for j in 1:npeople]
    
    #   Histogram of degrees
    println("mean number of neighbours:", mean(degrees))

    # histogram(degrees)
    # savefig("degree_distribution.png")
    # println("saving adjacency file")
    # for j in 1:npeople
    #     adjacency[j] = Int64.(adjacency[j])
    # end

    # open("adj.json", "w") do fp
    #     JSON.print(fp, adjacency)
    # end
    # open("positions.json","w") do fp
    #     JSON.print(fp, positions)
    # end
    return adjacency
end


function assign_neighbours!(agents::Vector{Agent})
    #   computes the adjacency matrix, representing the relationships between agents
    adjacency = run(agents)

    #   For each agent, it creates a new list of its neighbours - iterate over indicies j in agents adjacency list
    for (i, agent) in enumerate(agents)
        
        #   List of neighbours assigned to the neighbours attribute of agent
        agent.neighbours = [agents[j] for j in adjacency[i]]
    end
end


end