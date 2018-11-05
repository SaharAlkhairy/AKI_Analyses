function [] = analyseNRI_GEE()
% Description:
% computes the median NRI and pvalue across the different random
% number seeds;
% Inputs:
% GEE_NRI_IDI/NRI_IDI_GEE_models_Table_All_FludBal_noInOut_NoWHG_rng_{*}

TimeThresholdRange = [2:2:24];
NRI_IDI_GEE_models_AllRandomSeedsTable = table();
NRI_IDI_GEE_models_AllRandomSeedsTable.TimeThreshold  = TimeThresholdRange';

randomSeedNumList = 1:10; % used 1,2,3  % if just want to

for randomSeedNum = randomSeedNumList
    rng(randomSeedNum) % set seed for randomizer  in selecting the training and test sets in GEE models
    
    NRI_IDI_GEE_models_TableRaw = readtable("GEE_NRI_IDI\NRI_IDI_GEE_models_Table_All_FludBal_noInOut_NoWHG_rng_" + randomSeedNum + ".xlsx");
    
    % select columns TimeThreshold, NRI, NRI_pval: add suffix to NRI based
    % in random number seed
    
    NRI_IDI_GEE_models_Table = NRI_IDI_GEE_models_TableRaw(:,{'TimeThreshold', 'NRI', 'NRI_pval'}); 
    % change NRI and NRI_pval column names to _{randomSeedNum} which would
    % resolve the issue of join/ innerjoin changing the column names that
    % are the same to the name of the variable
    
    NRI_IDI_GEE_models_Table.Properties.VariableNames = {'TimeThreshold',  strcat('NRI_',num2str(randomSeedNum)), strcat('NRI_pval_' , num2str(randomSeedNum))};
    
    NRI_IDI_GEE_models_AllRandomSeedsTable = join(NRI_IDI_GEE_models_AllRandomSeedsTable, NRI_IDI_GEE_models_Table, 'Keys', ['TimeThreshold']);
   
    
end


   NRI_ColInds =  find(cellfun('length',regexpi(NRI_IDI_GEE_models_AllRandomSeedsTable.Properties.VariableNames, "NRI_[1-9*]")));
   NRI_pvalColInds = find(cellfun('length',regexpi(NRI_IDI_GEE_models_AllRandomSeedsTable.Properties.VariableNames, "NRI_pval_[1-9*]")));
   
  
   NRI_IDI_GEE_models_AllRandomSeedsTable.medianNRI = nanmedian(table2array(NRI_IDI_GEE_models_AllRandomSeedsTable(:,NRI_ColInds)), 2);
   NRI_IDI_GEE_models_AllRandomSeedsTable.medianNRI_pval = nanmedian(table2array(NRI_IDI_GEE_models_AllRandomSeedsTable(:,NRI_pvalColInds)), 2);
    
   writetable( NRI_IDI_GEE_models_AllRandomSeedsTable, "GEE_NRI_IDI\NRI_IDI_GEE_models_AllRandomSeedsTable.xlsx");
end

