%% Description:
% % All of this is from the selected cohort of 9621 samples 

% generates the figure 1 in manuscipt that plots the sensitivity and
% specificity for each of the time and volume combinations - curves
% generates the performances for each time and volume combination

% To add:
% compute number of samples (ICU IDs) that met the various time and volume thresholds
% compute NRI, p-value for NRI

% NewData matrix generated from UrineOutpurAnalysis_150604 
% NewData = [ICUSTAY_ID_COL GoldStandardColumnMat AKIUOCurrentdef];
% num columns of AKIUOCurrentdef = length(TimeThresholdRange) *
% length(AvgVolumeThresholdRange)


%% Note: changed the downloaded NetReclassificationImprovement and
% IntegratedDisciminationImprovement to include only indices in which
% neither the original nor the new prediction have NAs


%% Suggestions for plotting
% How about reducing eavch of the colored lines into a single point, that is the closest to the J-Point.
% 
% You might play around with adding gridlines.
% 
% You should follow – making pretty plots tutorial online to adjust the text style and size. 
% 
% You might want to make the dashed middle line thicker and a color which does not conflict with one of the other lines.
% 
% Not v, volume
% 
%% Inputs:
% TimeThresholdRange
% AvgVolumeThresholdRange
% NewData


%clear all
% clc

%% Plot properties
global width height alw fsz lw msz;
width = 3;% Width in inches
height = 3;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 18;%11;      % Fontsize
lw = 2; %1.5;      % LineWidth
msz = 8;       % MarkerSize


%% Loading 
load('TimeThresholdRange');
load('AvgVolumeThresholdRange');
load('NewData.mat');
addpath('..\NetReclassicationImprovement')
numIndices = length(TimeThresholdRange) * length(AvgVolumeThresholdRange);

%% calculate the specificty and sensitivity for each TimeThreshold and AvgVolumeThreshold combo

colnames_PerfTable =  {'linear_index','TimeThreshold','AvgVolumeThreshold','specificity','sensitivity', 'TN', 'FP', 'TP', 'FN', 'Sum_TN_FP_TP_FN', 'numMetThresh', 'dist100PercSensSpec', 'NRI', 'NRI_pval' };
PerformanceArray = NaN*ones(numIndices,length(colnames_PerfTable));% linear_index TimeThreshold AvgVolumeThreshold specificity sensitivity TN FP TP FN total
Performance = array2table(PerformanceArray,  'VariableNames', colnames_PerfTable);
   
% for NRI calculation  - index of 6, 0.5 prediction in NewData matrix
Ind_pred__6_0_5 =  2 + sub2ind([length(AvgVolumeThresholdRange),length(TimeThresholdRange)], find(AvgVolumeThresholdRange == 0.5), find(TimeThresholdRange == 6));
pred__6_0_5 =  NewData(:,Ind_pred__6_0_5);

GoldStandardColumnMat = NewData(:,2);

colSub =0;
for TimeThreshold = TimeThresholdRange
    colSub = colSub +1; 
    
    timeThreshInd = colSub;
    rowSub =0;
    for AvgVolumeThreshold = AvgVolumeThresholdRange
        rowSub=rowSub+1;
        volThreshInd = rowSub;
        % get Linear index / get column number
        linearInd = sub2ind([length(AvgVolumeThresholdRange),length(TimeThresholdRange)], rowSub, colSub);
        colNum = linearInd+2;
        
        % specify ICU
        %         %TFMatchesICU = (strcmp('MICU',ICUColumn) | strcmp('FICU',ICUColumn));
        %         TFMatchesICU = (strcmp('MICU',ICUColumn) | strcmp('FICU',ICUColumn));%strcmp('CCU',ICUColumn) ;
        %         RowsMatchingICU = find(TFMatchesICU==1);
        %         [CommonRows,ia,ib] = intersect(RowsMatchingICU,GoodRowsindices);
        
        % All
        
        % calculate confusion matrix
        [CM,order] = confusionmat(NewData(:,2),NewData(:,colNum));
        % calculate specificity = TN/(TN + FP)
        TN = CM(1,1);
        FP = CM(1,2);
        specificity = TN/(TN + FP);
        
        % calculate sensitivity = TP/(TP + FN)
        TP = CM(2,2);
        FN = CM(2,1);
        
        sensitivity = TP/(TP + FN);
        
        % TimeThreshold AvgVolumeThreshold specificity sensitivity TN FP TP FN
        Performance.linear_index(linearInd) = linearInd;
        Performance.TimeThreshold(linearInd) = TimeThreshold ;
        Performance.AvgVolumeThreshold(linearInd) = AvgVolumeThreshold ;
        Performance.specificity(linearInd) = specificity ;
        Performance.sensitivity(linearInd) = sensitivity ;
        Performance.TN(linearInd) = TN ;
        Performance.FP(linearInd) = FP ;
        Performance.TP(linearInd) = TP ;
        Performance.FN(linearInd) = FN ;
        Performance.Sum_TN_FP_TP_FN(linearInd) = sum([TN,FP,TP,FN]);
        Performance.numMetThresh(linearInd)= sum(Performance.TP(linearInd) + Performance.FP(linearInd) );
        Performance.dist100PercSensSpec(linearInd) =  sqrt((1 - Performance.specificity(linearInd))^2 + (1 -   Performance.sensitivity(linearInd))^2);
        
        % computing NRI values, and its pvalues
        timeVolThreshInd = 2 + sub2ind([length(AvgVolumeThresholdRange),length(TimeThresholdRange)], volThreshInd, timeThreshInd); % shifted because ICU ID  and CR gold standard column
        pred_time_volNew =  NewData(:, timeVolThreshInd);
        
        [ NRI , NRI_pval ] = NetReclassificationImprovement(pred__6_0_5, pred_time_volNew,GoldStandardColumnMat);
        % [ NRI , pval ] = NetReclassificationImprovement( pred_old , pred_new , outcome )
        Performance.NRI(linearInd) = NRI;
        Performance.NRI_pval(linearInd) = NRI_pval;
        
        Performance.sigNRI(linearInd) = Performance.NRI_pval(linearInd) <= 0.05;
    end
end

writetable(Performance, "AllTimeVolThresh_SensSpec_Dist_NRI.xlsx");
% save the performance table
sortPerfDist = sortrows(Performance, 'dist100PercSensSpec');


minTimes =   [2:2:24]; %[6,sortPerfDist.TimeThreshold(1:3)'];
minVols = [0 0.3 0.5 0.6 0.9];; % [0.5, sortPerfDist.AvgVolumeThreshold(1:3)'];
% find the best perfoming Vth and Tth combination

%[ mins]  = findMinDistPoint( Performance, 30);  % can't be used any more because performance is a table and not a matrix now.
%save('PerformanceAllCombo.mat','Performance');

% suggested improvement - see if can list time thresholds interested in and
% use contains? 

timesEvery2Indices =  find(ismember(Performance.TimeThreshold(:),minTimes ));



VolThr = minVols; %[0 0.3 0.5 0.6 0.9]; %[0 0.3 0.4 0.5 0.7 0.9];
colorShape = {'yx-','kx-','gx-','rx-','bx-'};

figure(2);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

for j =  1: length(VolThr)
    indicesVol = find(Performance.AvgVolumeThreshold(:) <VolThr(j) + 0.01 &Performance.AvgVolumeThreshold(:)> VolThr(j) - 0.01);
    indxs = intersect(timesEvery2Indices,indicesVol);
    
    
    plot(1-Performance.specificity(indxs),Performance.sensitivity(indxs),colorShape{j},'LineWidth',lw,'MarkerSize',msz);%,Performance(:,1)) %scatter
    
    hold on
end

%AvgVolumeThresholdRange = [0: 0.1: 1];%rowSub
%TimeThresholdRange = [2,4,6,8,10,12,14,16,18,20,22,24];% colSub



 % todo; if needed: make increase font size,make labels as (T= .., V = ..),
 % have automatic arrow.
ind6_pt5 = sub2ind([length(AvgVolumeThresholdRange),length(TimeThresholdRange)],6,3)
ind24_pt5 = sub2ind([length(AvgVolumeThresholdRange),length(TimeThresholdRange)],6,12)
ind2_pt5 = sub2ind([length(AvgVolumeThresholdRange),length(TimeThresholdRange)],6,1)
s = blanks(3)';
labels= [s,num2str(Performance.TimeThreshold([ind6_pt5,ind24_pt5,ind2_pt5])),s,num2str( Performance.AvgVolumeThreshold([ind6_pt5,ind24_pt5,ind2_pt5]))];

H = text(1- Performance.specificity([ind6_pt5,ind24_pt5,ind2_pt5]),Performance.sensitivity([ind6_pt5,ind24_pt5,ind2_pt5]),labels); %Performance(:,1)

set(H,'fontsize', 14)

xlim([0,1])
ylim([0,1])
plot([0,1],[0,1],'LineStyle','--','Color',[0.5 0.5 0.5],'LineWidth',lw)
xlabel('1 - Specificity')
ylabel('Sensitivity')
legend({'Volume = 0','Volume = 0.3','Volume = 0.5','Volume = 0.6','Volume = 0.9','Linear'});

grid on
%title('ROC plots')
legend('Location', 'SouthEast');

print('ROC plots','-dpng','-r300');


% Labeling all points
% s = blanks(numIndices)';
% labels= [s,num2str(Performance(:,2)),s,num2str( Performance(:,3))];
% plot(1-Performance(:,4),Performance(:,5),'x');%,Performance(:,1)) %scatter
% text(1- Performance(:,4),Performance(:,5),labels); %Performance(:,1)
% xlabel('1 - specificity')
% ylabel('sensitivity')
%title('MICU + FICU');
%zlabel('Linear Index')


% Specificity
% Specificity=true negatives/(true negative + false positives)
% If a person does not have the disease how often will the test be negative (true negative rate)?
%In other terms, if the test result for a highly specific test is positive you can be nearly certain that they actually have the disease.

%Sensitivity
% If a person has a disease, how often will the test be positive (true positive rate)?
% % Put another way, if the test is highly sensitive and the test result is negative you can be nearly certain that they don’t have disease.
% A Sensitive test helps rule out disease (when the result is negative). Sensitivity rule out or "Snout"
% Sensitivity= true positives/(true positive + false negative)


