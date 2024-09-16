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

"""
    Config

Configuration object
"""
@option struct Config
    num_datasets::Int = 1
    expression_method::String = "mean"
    file::FileConfig = FileConfig()
    format = FormatConfig()
    dataset::DatasetConfig = DatasetConfig()
    ontology::OntologyConfig = OntologyConfig()
    prefix::PrefixConfig = PrefixConfig()
    column::ColumnConfig = ColumnConfig()
end

"""
    path(config::Config)
"""
path(config::Config) = config.file.data_path

"""
    output_dir(config::Config)
"""
output_dir(config::Config) = config.file.output_dir

"""
    output_file(config::Config)
"""
output_file(config::Config) = joinpath(output_dir(config),
                                       config.dataset.base_eset_id *
                                       ".synthetic.jld2")

"""
    load_config(file_path::String)::Config

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
function load_config(file_path::String)::Config
    if isempty(file_path)
        @debug "No configuration file provided. Using default configuration."
        return Config()
    end

    try
        config_dict = TOML.parsefile(file_path)

        return from_dict(Config, config_dict)
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