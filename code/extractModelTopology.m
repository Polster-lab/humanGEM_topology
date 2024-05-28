%extractModelTopology
%
% Creates different binary matrices connecting the different molecular
% components in humanGEM  as well as metabolic subsytems
% (genes, metabolites and reactions). Results are stored as tab separated
% files in the results subfolder.
%
% Usage: set this as WD and type extractModelTopology in the matlab command 
% window.

%clone latest version of humanGEM
git('clone --depth=1 https://github.com/SysBioChalmers/Human-GEM.git')
load("Human-GEM/model/Human-GEM.mat")
mkdir('../results')
modelVer = ihuman.version;
dateStr  = datetime('today');
%get rxnGene matrix, model genes and stoichiometric matrix
[grRules,rxnGeneMat] = standardizeGrRules(ihuman);
genes = ihuman.genes;
short = ihuman.geneShortNames;
S = ihuman.S;
%get a list of ll unique metabolic subsystems in the model
subSystems = [];
for i=1:numel(ihuman.rxns)
    subS = ihuman.subSystems{i};
    subSystems = [subSystems;subS];
end
subSystemsAll = subSystems;
subSystems = unique((subSystems));
SS = strrep(subSystems,' ','_');
SS = strrep(SS,'-','_');
SSs = {};
for i=1:length(subSystems)
    SSs{i} = ['SS_' num2str(i)];
end
%write a file with the rxnGeneMat
c = array2table(S);
t = table(ihuman.mets);
t = [t,c];
t.Properties.VariableNames = [{'mets'},ihuman.rxns'];
writetable(t,'../results/S_matrix.txt','delimiter','\t','QuoteStrings',false)
%write a file with the rxnGeneMat
c = array2table(rxnGeneMat);
t = table(ihuman.rxns);
t = [t,c];
t.Properties.VariableNames = [{'rxns'},ihuman.genes'];
writetable(t,'../results/rxnGeneMatrix.txt','delimiter','\t','QuoteStrings',false)
%obtain a binary matrix indicating the presence of metabolites in reactions
%catayzed by each of the model's genes [RGM x logical(S)']
GeneMetsMatrix = full(rxnGeneMat'*logical(S'));
c = array2table(GeneMetsMatrix);
t = table(ihuman.genes);
t = [t,c];
t.Properties.VariableNames = [{'genes'},ihuman.mets'];
writetable(t,'../results/geneMetMatrix.txt','delimiter','\t','QuoteStrings',false)
%generate a matrix indicating the metabolic subsystem for reactions
%catalysed by each gene's products
rxnSubsystMat = zeros(length(ihuman.rxns),length(subSystems));
geneSubsystMat = zeros(length(ihuman.genes),length(subSystems));
metSubSystMat = zeros(length(ihuman.mets),length(subSystems));
for i=1:numel(subSystems)
    x = find(contains(subSystemsAll,subSystems{i}));
    for j=1:length(x)
        genes = find(rxnGeneMat(x(j),:));
        mets  = find(S(:,i));
        geneSubsystMat(genes,i) = 1;
        metSubSystMat(mets,i) = 1;
    end
end
c = array2table(geneSubsystMat);
c.Properties.VariableNames = SSs;
t = table(ihuman.genes);
t = [t,c];
t.Properties.VariableNames(1) = {'genes'};
writetable(t,'../results/geneSubSystemMatrix.txt','delimiter','\t','QuoteStrings',false)
%
c = array2table(metSubSystMat);
c.Properties.VariableNames = SSs;
t = table(ihuman.mets);
t = [t,c];
t.Properties.VariableNames(1) = {'mets'};
writetable(t,'../results/metSubSystemMatrix.txt','delimiter','\t','QuoteStrings',false)
%get lists with rxns ids names and formulas and grRules, also another one
%with met unique IDs and their names as well compartmentalization
formulas = constructEquations(ihuman);
t = table(ihuman.rxns,ihuman.rxnNames,formulas,grRules);
t.Properties.VariableNames = {'rxns' 'rxnNames' 'formulas' 'grRules'};
writetable(t,'../results/ihuman_rxns.txt','delimiter','\t','QuoteStrings',false)
%
t = table(ihuman.mets,ihuman.metNames,ihuman.compNames(ihuman.metComps),ihuman.metFormulas);
t.Properties.VariableNames = {'mets' 'metNames' 'compartment' 'metFormulas'};
writetable(t,'../results/ihuman_mets.txt','delimiter','\t','QuoteStrings',false)
%print a file with date and model version
fileID = fopen('../results/model_version.txt','w');
formatSpec = 'ihuman version: %s\ndate: %s';
fprintf(fileID,formatSpec,modelVer,dateStr)
fclose(fileID);
rmdir('Human-GEM','s')