using SpaceGroups
using ChargeFlipPhaser

function load_data(file_path::String)
    # TODO: Load the data from JSON file and create DiffractionData object as
    # well as the Vector{Complex{Float64}} for the structure factors.

    # When recrationg the space groups, make a sanity check - when all group elements are
    # passed as generators, the order of the group should be equal to the number of elements passed
    # (no other elements should appear in the group).

    # Just to check the dependencies, erase these lines later:
    e1 = @SGE([0 1; -1 0]);
    e2 = @SGE([1//1, 1//1]);
    println(e1∘e2)

end

filepath=joinpath(@__DIR__, "data", "CdYb.json")
load_data(filepath)