
module SyntheticExpressionMixtures

using Comonicon
using Configurations
using DataFrames
using Distributions
using ExpressionData
using UUIDs

using TOML

include("./config.jl")
include("./utils.jl")
include("./functions.jl")

"""
    sem

Generate synthetic i.e. *in silico* expression mixtures from a given base expression set.

# Args
- `eset_file`: The path to the base expression set file.

# Options
- `-c, --config_file`: The path to the configuration file. Default is `config.toml`.
"""
@main function sem(eset_file::String; config_file::String="")
    # Note that if the file is not found, the default configuration will be used.
    # - this is a feature of the Configurations.jl package 
    config::Config = load_config(config_file)

    # Load the input expression data
    input::ExpressionSet = load_input(eset_file, config)

    # Generate the synthetic expression mixtures
    synthetic_eset, descriptor = generate_synthethic_expression_mixtures(input,
                                                                         config)

    return save_results(synthetic_eset, descriptor, config)
end

end
