
module SyntheticExpressionMixtures

using Configurations
using DataFrames
using Distributions
using ExpressionData
using UUIDs

include("./config.jl")

include("./utils.jl")
export generate_proportions, calculate_gene_expressions, create_proportion_label,
       generate_descriptor_file

include("./functions.jl")
export generate_synthetic_expression_values

end
