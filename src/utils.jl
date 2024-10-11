"""
    generate_proportions(numclasses::Int; samples::Int=100, decimals::Int=2, id_prefix="proportion_")

Generate random proportions for each class. The proportions are generated as a DataFrame with
`samples` rows and `numclasses` columns. The column names are generated using the `id_prefix` and
the column index.

The proportions are generated as random numbers between 0 and 1, and then scaled so that their sum
equals 100. The proportions are rounded to the specified number of decimal places.
"""
function generate_proportions(N::Int;
                              samples::Int=100,
                              decimals::Int=2,
                              id_prefix="proportion_")
    # Create a Dirichlet distribution with N classes
    dirichlet = Dirichlet(repeat([1.0], N))

    # Sample N proportions from the Dirichlet distribution
    proportions = rand(dirichlet, samples)

    # Scale the proportions so that they sum to 100
    columns = [Symbol("$(id_prefix)$i") for i in 1:samples]

    return DataFrame(proportions, columns)
end

"""
    calculate_gene_expressions(cell_data::Matrix, method::Symbol)

Calculate the gene expression values for a given cell type based on the method provided.
"""
function calculate_gene_expressions(cell_data::Matrix, method::String)
    if method == "sum"
        return sum(cell_data; dims=2)
    elseif method == "mean"
        return mean(cell_data; dims=2)
    else
        throw(ArgumentError("Unknown expression method: $method"))
    end
end

"""
    create_proportion_label(cell_types::Vector{String}, proportions::Vector{Float64};
                           cell_separator=";", value_separator="=", value_unit="%")

Create a label string representing the proportions of different cell types. The label is
formatted as `cell_type1<<value_seperator>>value1<<value_unit>><<cell_separator>>cell_type2<<value_seperator>>..."

By default, the format is `cell_type1=value1%;cell_type2=value2%...`


"""
function create_proportion_label(cell_types::Vector{String},
                                 proportions::Vector{Float64};
                                 cell_separator=";",
                                 value_separator="=",
                                 value_unit="%")
    return join([create_proportion_label(cell_types[i], proportions[i];
                                         value_separator=value_separator,
                                         value_unit=value_unit)
                 for i in eachindex(cell_types)],
                cell_separator)
end

function create_proportion_label(cell_type::String,
                                 proportion::Float64;
                                 value_separator="=",
                                 value_unit="%")
    return "$cell_type$value_separator$proportion$value_unit"
end

"""
    generate_descriptor_file(config::SYDConfig)
Generate a descriptor file for the synthetic datasets. The descriptor file is a YAML file	
containing the metadata of the synthetic datasets. The metadata includes the dataset ID,
title, and the series of synthetic datasets. The series includes the ID, platform, and type of each dataset.
"""
function generate_descriptor_file(config::SYDConfig)

    # Get the configuration values
    (; base) = config.prefix
    (; samples, name_template, base_eset_id) = config.dataset
    num_datasets = config.num_datasets

    dataset_id = "$(base)_$(uuid4())"
    dataset_name = replace(name_template, "N" => string(samples),
                           "BASE_ESET_ID" => base_eset_id)

    # There is only one series in this case                           
    series = [Dict("id" => dataset_id,
                   "platform" => "synthetic",
                   "type" => "expression")]

    descriptor = Dict(dataset_id => Dict("title" => dataset_name,
                                         "id" => dataset_id,
                                         "series" => series))

    return descriptor
end
