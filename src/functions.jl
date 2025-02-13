"""
    identify_cell_types(eset::ExpressionSet, config::SYDConfig)

Identify the cell types in the expression set based on the configuration provided.
The function assumes that the expression set has a column with the cell type information for each
sample. The cell type column is specified in the config.
"""
function identify_cell_types(eset::ExpressionSet,
                             config::SYDConfig)::Dict{String,Vector{Int}}
    # Get the cell type column
    (; cell_type) = config.column

    # Get the cell types
    cell_types = phenotype_data(eset)[!, cell_type]

    # Find the indices of the samples for each cell type
    # This builds a dictionary with the cell type as the key and the indices of the samples in the
    # original expression set as the value
    cell_type_indices = Dict{String,Vector{Int}}()
    for cell_type in unique(cell_types)
        cell_type_indices[cell_type] = findall(==(cell_type), cell_types)
    end

    return cell_type_indices
end

"""
    aggregate_expression_values(eset::ExpressionSet, cell_types::Dict{String,Vector{Int}}, config::SYDConfig)
"""
function aggregate_expression_values(eset::ExpressionSet,
                                     cell_types::Dict{String,Vector{Int}},
                                     config::SYDConfig)::ExpressionSet
    # Get the expression method
    (; expression_method) = config

    # Get the expression values
    gxdata = expression_values(Matrix, eset)

    # To get the new expression, lets reduce the original expression based on the deried methods
    # and the target proportions
    reduced_gxdata = [calculate_gene_expressions(gxdata[:, cell_types[cell_type]],
                                                 expression_method)
                      for cell_type in keys(cell_types)]
    reduced_gxdata = hcat(reduced_gxdata...)

    # Return the reduced expression set
    return ExpressionSet(reduced_gxdata,
                         phenotype_data(eset),
                         feature_data(eset),
                         experiment_data(eset),
                         annotation(eset))
end

"""
    generate_synthetic_expression_values(eset::ExpressionSet, cell_types::AbstractVector{String},
                                         proportions::DataFrame, config::SYDConfig)

Generate synthetic expression values for the given expression set based on the cell types, proportions,
and configuration provided.

The function will generate the synthetic expression values for each sample based on the proportions
and the aggregated expression values for each cell type.

The synthetic expression values are calculated as follows:
1. For each sample, iterate over the cell types and their corresponding proportions.
2. For each cell type, calculate the target expression value by multiplying the aggregated expression
   value for the cell type by the target proportion.
3. Sum the target expression values for all cell types to get the synthetic expression value for the
    sample.

The function returns a matrix of synthetic expression values with the same number of rows as the
original expression set and the number of columns equal to the number of samples.
"""
function generate_synthetic_expression_values(eset::ExpressionSet,
                                              cell_types::Dict{String,Vector{Int}},
                                              proportions::DataFrame,
                                              config::SYDConfig)::Matrix
    (; noise) = config
    noise_gen = noise.method == "normal" ? Normal(noise.mean, noise.std) :
                Uniform(noise.min, noise.max)

    # Get the expression values
    gxdata = expression_values(Matrix, eset)

    # The new_expression will be the same as gxdata but with fewer samples
    new_expression = zeros(size(gxdata, 1), size(proportions, 2))

    # To get the new expression, lets reduce the original expression based on the deried methods
    # and the target proportions
    for (sample_index, sample) in enumerate(names(proportions))
        for (cell_type_index, cell_type) in enumerate(keys(cell_types))
            target_proportion = proportions[!, sample][cell_type_index]

            cell_type_indices = cell_types[cell_type]

            if length(cell_type_indices) > 1
                gene_exprs = gxdata[:, cell_type_index]
                gx_ = vec(gene_exprs) .* target_proportion
                if noise.noise
                    new_expression[:, sample_index] += gx_ .+ rand(noise_gen, size(gx_))
                else
                    new_expression[:, sample_index] += gx_
                end
            else
                if noise.noise
                    new_expression[:, sample_index] += gxdata[:,
                                                              first(cell_type_indices)] .*
                                                       target_proportion .+
                                                       rand(noise_gen, size(gxdata, 1))
                else
                    new_expression[:, sample_index] += gxdata[:,
                                                              first(cell_type_indices)] .*
                                                       target_proportion
                end
            end
        end
    end

    return new_expression
end

function generate_synthetic_mixture_pdata(base_eset::ExpressionSet,
                                          proportions::DataFrame,
                                          config::SYDConfig)::DataFrame
    # Get the cell type column
    (; cell_type) = config.column

    # Get the cell types
    cell_types = unique(phenotype_data(base_eset)[!, cell_type])

    # Create the proporiton label 
    proportions_columns = Array{String}(undef, size(proportions, 2))
    for (sample_index, sample) in enumerate(names(proportions))
        proportions_columns[sample_index] = create_proportion_label(cell_types,
                                                                    proportions[!, sample])
    end

    return DataFrame(; sample_names=names(proportions),
                     cell_type=proportions_columns)
end

"""
    generate_synthetic_expression_mixtures(base_eset::ExpressionSet,
                                           config::SYDConfig)

Generate synthetic expression mixtures from a given base expression set. The function will generate
the synthetic expression mixtures based on the configuration provided.
"""
function generate_synthetic_expression_mixtures(base_eset::ExpressionSet,
                                                config::SYDConfig)

    # Step 1: Sample the synthetic proportions corresponding to each cell types
    #       This function will generate the synthetic proportions for each cell type
    #       The function will return a DataFrame with the synthetic proportions
    # -- relevant configs:
    #       - samples
    num_classes = length(unique(phenotype_data(base_eset)[!, config.column.cell_type]))
    proportions = generate_proportions(num_classes; samples=config.dataset.samples)

    # Step 2: Identify the cell types in the expression set
    #       This step assumes that the expression set has a column with the cell type information
    #       for each sample, if this is not the case, the program will exit with an error
    #       The cell type column is specified in the configuration file
    #       On return, the function will return a dictionary with the cell types as keys and the
    #       indices of the samples for each cell type as values
    # -- relevant configs:
    #       - cell_type_column
    cell_types::Dict{String,Vector{Int}} = identify_cell_types(base_eset, config)

    # Step 3: Aggregate the expression values for each cell type using the specified method
    #       The method is specified in the configuration file
    #       This function will return a new expression set with the aggregated expression values
    #       for each cell type
    # -- relevant configs:
    #       - expression_method
    base_eset_aggr::ExpressionSet = aggregate_expression_values(base_eset, cell_types,
                                                                config)

    # Step 4: Generate the synthetic expression values for each sample
    #       This function will generate the synthetic expression values for each sample based on the
    #       proportions and the aggregated expression values for each cell type
    syn_gxdata::Matrix = generate_synthetic_expression_values(base_eset_aggr,
                                                              cell_types,
                                                              proportions, config)

    # Step 5: Generate the synthetic expression set phenotype data 
    #       This function will generate the synthetic expression set (structs) based on the synthetic
    #       expression values and the original expression set
    sym_pdata::DataFrame = generate_synthetic_mixture_pdata(base_eset, proportions, config)

    # Create the synthetic expression set
    synthetic_eset = ExpressionSet(syn_gxdata,
                                   sym_pdata,
                                   feature_data(base_eset),
                                   experiment_data(base_eset),
                                   annotation(base_eset))

    return synthetic_eset
end

"""
    generate_synthetic_expression_mixtures(config::SYDConfig)

Generate a random synthetic expression set based on the configuration provided. This function will
generate a synthetic expression set with random expression values and proportions based on the
configuration provided.
"""
function generate_synthetic_expression_mixtures(config::SYDConfig)
    # Step 1: Generate the synthetic expression set
    #       This function will generate a synthetic expression set with random expression values
    #       and proportions based on the configuration provided
    # -- relevant configs:
    #       - num_genes
    #       - num_samples
    #       - num_classes
    #       - samples
    base_eset = generate_synthetic_expression_set(config)

    # Step 2: Generate the synthetic expression mixtures
    #       This function will generate the synthetic expression mixtures based on the configuration
    #       provided
    # -- relevant configs:
    #       - cell_type_column
    #       - expression_method
    #       - noise
    synthetic_eset = generate_synthetic_expression_mixtures(base_eset, config)

    return synthetic_eset
end

"""
    generate_synthetic_expression_set(config::SYDConfig)

Generate a synthetic expression set with random expression values. 
"""
function generate_synthetic_expression_set(n_samples, n_genes;)
    # Generate the expression values
    gxdata = randn(n_genes, n_samples)

    # Generate the phenotype data
    sample_names = ["sample_$i" for i in 1:n_samples]
    cell_types = ["cell_type_$i" for i in 1:n_samples]
    pdata = DataFrame(; sample_names=sample_names, cell_type=cell_types)

    # Generate the feature data
    feature_names = ["gene_$i" for i in 1:n_genes]
    fdata = DataFrame(; feature_names=feature_names)

    # Generate the experiment data
    edata = DataFrame(; experiment_data=Dict{String,Any}())

    # Generate the annotation
    annotation = DataFrame(; annotation=Dict{String,Any}())

    return ExpressionSet(gxdata, pdata, fdata, edata, annotation)
end