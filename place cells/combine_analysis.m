function analysis = combine_analysis(analysis1, analysis2)
% Combine two analysis structures
% Useful for combining multiplane data

if isempty(analysis1) && isempty(analysis2)
    error('WTF?!');
elseif isempty(analysis1)
    analysis = analysis2;
    return
elseif isempty(analysis2)
    analysis = analysis1;
    return
end

analysis = analysis1;
analysis.psth = cat(3, analysis1.psth, analysis2.psth);
analysis.raw_psth = cat(3, analysis1.raw_psth, analysis2.raw_psth);
analysis.raw_stack = cat(2, analysis1.raw_stack, analysis2.raw_stack);
analysis.Pi = cat(2, analysis1.Pi, analysis2.Pi);
analysis.stack = cat(2, analysis1.stack, analysis2.stack);
analysis.SI = cat(2, analysis1.SI, analysis2.SI);
analysis.pval = cat(2, analysis1.pval, analysis2.pval);
analysis.sparsity = cat(2, analysis1.sparsity, analysis2.sparsity);
analysis.width = cat(2, analysis1.width, analysis2.width);
analysis.original_deconv = cat(2, analysis1.original_deconv, analysis2.original_deconv);
analysis.silent = cat(2, analysis1.silent, analysis2.silent);
analysis.rick_rejects = cat(2, analysis1.rick_rejects, analysis2.rick_rejects);
analysis.SI_marge = cat(2, analysis1.SI_marge, analysis2.SI_marge);
analysis.pc_list = cat(2, analysis1.pc_list, analysis2.pc_list + length(analysis1.psth));

try analysis = rmfield(analysis, 'deconv'); catch; end
try analysis = rmfield(analysis, 'shuff_stack'); catch; end