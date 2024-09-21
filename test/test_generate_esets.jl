
@testsnippet MockEset begin
    using ExpressionData

    cell_types = ["B cells", "T cells", "NK cells", "Monocytes", "Dendritic cells"]
    eset = rand(ExpressionSet, 5, 10)
end

@testsnippet MockConfig begin
    config = Config()
end

@testitem "generate expression mixtures " setup = [MockEset, MockConfig] begin
    generate_synthethic_expression_values(eset, config)
end