git('clone --depth=1 https://github.com/SysBioChalmers/Human-GEM.git')
load("Human-GEM/model/Human-GEM.mat")
[grRules,GeneRxnM] = standardizeGrRules(ihuman);
genes = ihuman.genes;
short = ihuman.geneShortNames;
S = ihuman.S;
GeneMetsMatrix = full(GeneRxnM'*logical(S'));
