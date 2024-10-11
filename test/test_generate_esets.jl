
@testsnippet MockEset begin
    using ExpressionData

    cell_types = ["B cells", "T cells", "NK cells", "Monocytes", "Dendritic cells"]
    eset = rand(ExpressionSet, 5, length(cell_types))
    eset.exprs .= rand(1:100, size(eset.exprs))
    eset.phenotype_data[!, :cell_types] = cell_types
end

@testsnippet MockConfig begin
    config = SyntheticExpressionMixtures.SYDConfig(;
                                                   column=SyntheticExpressionMixtures.ColumnConfig(;
                                                                                                   cell_type="cell_types"))
end

@testitem "generate expression mixtures " setup = [MockEset, MockConfig] begin
    @time sy_eset = generate_synthetic_expression_mixtures(eset, config)

    @test isa(eset, ExpressionSet)
end