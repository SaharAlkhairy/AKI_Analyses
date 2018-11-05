function [ Table1 ] = CreateTable1( data, AllRows, GoodRowsTF)
Table1 = cell(11,6);

for i = 1: 6
    if i==1 % all
        TFMatchesICU = GoodRowsTF;
        RowsMatchingICU = find(TFMatchesICU==1);
        numAll = length(RowsMatchingICU);
        
    end
    
    if i==2 % CSRU
        ICU = 'CSRU';
        ICUColumn = ConvertToCell(data,14,AllRows);
        TFMatchesICU = (strcmp(ICU,ICUColumn))& GoodRowsTF;
        RowsMatchingICU = find(TFMatchesICU==1);
        
        numCSRU = length(RowsMatchingICU);
    end
    
    if i==3 % SICU
        ICU = 'SICU';
        ICUColumn = ConvertToCell(data,14,AllRows);
        TFMatchesICU = (strcmp(ICU,ICUColumn))& GoodRowsTF;
        RowsMatchingICU = find(TFMatchesICU==1);
        numSICU = length(RowsMatchingICU);

    end
    
    if i==4 % CCU
        ICU = 'CCU';
        ICUColumn = ConvertToCell(data,14,AllRows);
        TFMatchesICU = (strcmp(ICU,ICUColumn))& GoodRowsTF;
        RowsMatchingICU = find(TFMatchesICU==1);
        numCCU = length(RowsMatchingICU);

    end
    
    if i==5 % MICU + FICU
        ICU1 = 'MICU';% 'MICU'+'FICU'
        ICU2 = 'FICU';
        ICUColumn = ConvertToCell(data,14,AllRows);
        TFMatchesICU = (strcmp(ICU1,ICUColumn) | strcmp(ICU2,ICUColumn))& GoodRowsTF;
        RowsMatchingICU = find(TFMatchesICU==1);
        numMICU = length(RowsMatchingICU);

    end
    
    if i == 6 % all population
        RowsMatchingICU = AllRows';
    end
    
    % number of patients
    NumPatients = size(RowsMatchingICU,1);
    if i==1
        TotalNumber = NumPatients;
    end
    % Number of patients
    Table1{1,i} = [num2str(NumPatients),' (',num2str((NumPatients*100)/TotalNumber),' %)'];
    
    % - Age
    MedianAge = nanmedian(cell2mat(ConvertToCell(data,27,RowsMatchingICU))); % after removing NaN:   62.2977
    Q = quantile(cell2mat(ConvertToCell(data,27,RowsMatchingICU)),[0.25, 0.75]);  %  17.8035
    
    Table1{2,i} = [num2str(MedianAge), ' (',num2str(Q(1)),', ',num2str(Q(2)),')'];
    
    % SAPSI - first
    MedianSAPSI_first = nanmedian(cell2mat(ConvertToCell(data,22,RowsMatchingICU))); %   13.1826
    Q = quantile(cell2mat(ConvertToCell(data,22,RowsMatchingICU)),[0.25, 0.75]);  %     5.3234
    Table1{3,i} = [num2str(MedianSAPSI_first), ' (',num2str(Q(1)),', ',num2str(Q(2)),')'];
    
    % SOFA - first
    MedianSOFA_first = nanmedian(cell2mat(ConvertToCell(data,21,RowsMatchingICU))); %
    Q = quantile(cell2mat(ConvertToCell(data,21,RowsMatchingICU)),[0.25, 0.75]);  %
    Table1{4,i} = [num2str(MedianSOFA_first), ' (',num2str(Q(1)),', ',num2str(Q(2)),')'];
    
    % Gender
    numMale = length(find(strcmp('M',ConvertToCell(data,15,RowsMatchingICU))==1)); %
    numFemale = length(find(strcmp('F',ConvertToCell(data,15,RowsMatchingICU))==1)); %
    percentMale = (numMale*100)/NumPatients;
    %Table1{5,i} = [num2str(percentMale), ' %'];
    Table1{5,i} = [num2str(numMale),' (',num2str(percentMale),' %)'];

    
    %   first hospital service
%     figure()
%     [idx,label] = grp2idx(ConvertToCell(data,14,RowsMatchingICU));
%     
%     hist(idx,unique(idx));
%     set(gca,'xTickLabel',label)
%     title('ICUSTAY FIRST SERVICE');
    
    % length of stay should only include deceased patients
    deseased = find(strcmp('Y',ConvertToCell(data,20,AllRows))==1);% return back to 'N' to have only survivals
    deseasedIntersect = intersect(deseased,RowsMatchingICU);
    
    alive = find(strcmp('N',ConvertToCell(data,20,AllRows))==1);
    aliveIntersect = intersect(alive,RowsMatchingICU);
    
    
    %     hospital length of stay - alive
    MedianHOSPITAL_LOS = nanmedian(cell2mat(ConvertToCell(data,17,aliveIntersect)))/(60*24); %   1.6700e+04
    Q = quantile(cell2mat(ConvertToCell(data,17,aliveIntersect)),[0.25,0.75])/(60*24);  %       2.0432e+04
    Table1{6,i} = [num2str(MedianHOSPITAL_LOS), ' (',num2str(Q(1)),', ',num2str(Q(2)),')'];
    
        %     hospital length of stay - deceased
    MedianHOSPITAL_LOS = nanmedian(cell2mat(ConvertToCell(data,17,deseasedIntersect)))/(60*24); %   1.6700e+04
    Q = quantile(cell2mat(ConvertToCell(data,17,deseasedIntersect)),[0.25,0.75])/(60*24);  %       2.0432e+04
    Table1{7,i} = [num2str(MedianHOSPITAL_LOS), ' (',num2str(Q(1)),', ',num2str(Q(2)),')'];
    
    %     - ICU length of stay -alive
    MedianICU_LOS = nanmedian(cell2mat(ConvertToCell(data,18,aliveIntersect)))/(60*24); %    7.3687e+03
    Q = quantile(cell2mat(ConvertToCell(data,18,aliveIntersect)),[0.25,0.75])/(60*24);  %        1.7593e+04
    Table1{8,i} = [num2str(MedianICU_LOS), ' (',num2str(Q(1)),', ',num2str(Q(2)),')'];
    
     %     - ICU length of stay - deceased
    MedianICU_LOS = nanmedian(cell2mat(ConvertToCell(data,18,deseasedIntersect)))/(60*24); %    7.3687e+03
    Q = quantile(cell2mat(ConvertToCell(data,18,deseasedIntersect)),[0.25,0.75])/(60*24);  %        1.7593e+04
    Table1{9,i} = [num2str(MedianICU_LOS), ' (',num2str(Q(1)),', ',num2str(Q(2)),')'];
    
    % Hospital survival rate
    numDeath = length(find(strcmp('Y',ConvertToCell(data,20,RowsMatchingICU))==1)); % 1348
    numAlive = length(find(strcmp('N',ConvertToCell(data,20,RowsMatchingICU))==1)); % 14401
    
    survival_rate = (numAlive * 100)/(numAlive + numDeath); %  0.9144
    
    Table1{10,i} = [num2str(survival_rate), ' %'];
    
    % First Weight
    MedianFirstICUWeight =nanmedian(cell2mat(ConvertToCell(data,68,RowsMatchingICU)));
    
    Q = quantile(cell2mat(ConvertToCell(data,68,RowsMatchingICU)),[0.25,0.75]);
    
    Table1{11,i} = [num2str(MedianFirstICUWeight), ' (',num2str(Q(1)),', ',num2str(Q(2)),')'];
    
    
end



width = 3;     % Width in inches
height = 3;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize


figure()

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
bar([numCSRU,numSICU,numCCU,numMICU]);
set(gca,'xTickLabel',{'CSRU','SICU','CCU','MICU'});
title('ICUSTAY FIRST SERVICE');
print('ICUSTAY FIRST SERVICE','-dpng','-r300');



end

