using Configurations

@option struct FileConfig
    data_path::String = "/data"
    output_dir::String = joinpath("/data", "eset.synthetic.jld2")
    datasets_descriptor_file::String = "datasets.yaml"
    generate_datasets_descriptor_file::Bool = true
end

@option struct FormatConfig
    decimals::Int = 2
    value_unit::String = "%"
    cell_separator::String = ";"
    value_separator::String = "="
end

@option struct OntologyConfig
    ontology_id::String = "CL"
    ontology_name::String = "cell type"
end

@option struct DatasetConfig
    name_template::String = "Synthetic dataset of N samples. Derived from BASE_ESET_ID" # Note the two variables available for substitution are N and BASE_ESET_ID
    base_eset_id::String = ""
    samples::Int = 1000
end

# Lets bundle the config properties related to nomeclature in a custom struct
@option struct PrefixConfig
    base::String = "SD"
    base_pad::Int = 3
    run::String = "SR"
    run_pad::Int = 4
    sample_id::String = "proportion_"
end

@option struct ColumnConfig
    cell_type::String = "cell type"
    proportion::String = "proportion"
    feature_id::String = "feature id"
end

# This config handles how the noise is added to the synthetic data
@option struct NoiseConfig
    noise::Bool = true
    method::String = "normal"
    mean::Float64 = 0.0
    std::Float64 = 0.1
    min::Float64 = 0.0
    max::Float64 = 1.0
end

"""
    SYDConfig

Configuration object
"""
@option struct SYDConfig
    num_datasets::Int = 1
    expression_method::String = "mean"
    file::FileConfig = FileConfig()
    format = FormatConfig()
    dataset::DatasetConfig = DatasetConfig()
    ontology::OntologyConfig = OntologyConfig()
    prefix::PrefixConfig = PrefixConfig()
    column::ColumnConfig = ColumnConfig()
    noise::NoiseConfig = NoiseConfig()
end

"""
    path(config::SYDConfig)
"""
path(config::SYDConfig) = config.file.data_path

"""
    output_dir(config::SYDConfig)
"""
output_dir(config::SYDConfig) = config.file.output_dir

"""
    output_file(config::SYDConfig)
"""
output_file(config::SYDConfig) = joinpath(output_dir(config),
                                          config.dataset.base_eset_id *
                                          ".synthetic.jld2")

create_config(conf::AbstractDict) = from_dict(SYDConfig, conf)

"""
    load_config(file_path::String)::SYDConfig

Load the configuration from a TOML file. It reads the TOML file and extracts the input strings,
output file, and verbose flag from the configuration. It returns a Config object with the extracted
values.

If the file is not found or the TOML format is invalid, it prints an error message and exits the program.

Example TOML configuration file:
```
output_file = "graph.png"
verbose = true
```
"""
function load_config(file_path::String)::SYDConfig
    if isempty(file_path)
        @debug "No configuration file provided. Using default configuration."
        return SYDConfig()
    end

    try
        config_dict = TOML.parsefile(file_path)

        return from_dict(SYDConfig, config_dict)
    catch e
        if isa(e, SystemError)
            @error "File not found: $file_path"
        elseif isa(e, TOML.ParserError)
            @error "Invalid TOML format in $file_path"
        else
            rethrow(e)
        end
        exit(1)
    end
end

const supported_eset_extensions = [".rds", ".RDS", ".jld2"]
"""
    load_input(eset_file::String, config::SYDConfig)

Load the input expression data from the file specified in the configuration. The function will
load the expression data from the file and return it as an `ExpressionSet`.
"""
function load_input(eset_file::String, config::SYDConfig)

    # Check if the file extension is valid
    if all(endswith.(eset_file, supported_eset_extensions) .== false)
        throw(ArgumentError("Invalid file extension: $eset_file. Supported extensions: $(join(supported_eset_extensions, ", "))"))
        exit(1)
    end

    # Load the expression data from the file
    eset = load_eset(eset_file)

    # Get the number of samples from the configuration
    samples = config.dataset.samples

    # If the number of samples is less than the number of samples in the expression set,
    # select a random subset of samples
    if samples < length(names(phenotype_data(eset)))
        indices = rand(1:size(eset, 2), samples)
        eset = eset[:, indices]
    end

    return eset
end

function save_results(eset::ExpressionSet, descriptor::Dict, config::SYDConfig)
    try
        # Writing the results to output dir
        eset_filename = output_file(config)

        if isfile(eset_filename)
            @warn "Overwriting existing expression set file: $eset_filename"
        end

        # Save the synthetic expression set to a file
        try
            save_eset(eset, eset_filename)
            @info "Saved results to: $eset_filename"
        catch e
            @error "Failed to save eset"
        end

        # Save the descriptor file
        descriptor_filepath = joinpath(output_dir(config), "descriptor.toml")

        # Create the file
        if isfile(descriptor_filepath)
            @warn "Overwriting existing descriptor file: $descriptor_filepath"
        else
            mkpath(dirname(descriptor_filepath))
        end

        # Populate file with data from dict using TOML
        open(descriptor_filepath, "w") do io
            return TOML.print(io, descriptor)
        end

        @info "Saved descriptor to: $descriptor_filepath"

    catch e
        if isa(e, SystemError)
            @error "Unable to access the file system: $e"
            exit(1)
        end
        @error "An error occured while saving: $e"
    end
end