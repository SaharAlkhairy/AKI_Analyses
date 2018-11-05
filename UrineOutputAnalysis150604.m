%% Toolboxes needed:
% deep learning toolbox
% Parallel Computing Toolbox.

projPath = "~\"  % edit this to path to this file
cd(projPath)  
%% Adding paths for methods and other - note: I edited these functions 
addpath('NetReclassificationImprovement\') 
addpath('IntegratedDiscriminationImprovement\')
addpath('bin\')

%% Using new AKICr
%% New changesTruncatingUOData
% removing the Subject ID only dealing with ICU ID: changed from 3+ lnIndx
% to 2+ linindex
% adding avg volumes per window
% Removing extraction of below threshold features
% Trying to vectorize window extraction

% _2
% trying to vectorize timeThresholds

% _3:
% Preallocated the GEE arrays with ~ 100 hours * 11050 size.

% _5: 
% ROC with CI and added LBM and with uptil 3 hours prior MAP instead of
% only 1. Beyond three hours might not be reliable. 

%cinner = join(a,b,'key','Key1','Type','inner', 'MergeKeys',true)
%% fixing the data set by only having non NaN and common ICUs
% indicesnonNaNVol = find(isnan(urinedata.VOLUME)==0);
% urinedataNonNaN = urinedata(indicesnonNaNVol,:);
%
% % remove NaN weight
% indicesnonNaNWeight = find(isnan(DataSetAllWithFirstCrMaxWeightNoRepeat.MAX_WEIGHT)==0);
%
% % getting only the common ICU IDs
%
% [c1,ia1,ib1] = intersect(unique(DataSetAllWithFirstCrMaxWeightNoRepeatCommon(AllRows,:).ICUSTAY_ID),unique(exportNew.ICUSTAY_ID));
% length(c1)
% DataSetAllWithFirstCrMaxWeightNoRepeatCommon = DataSetAllWithFirstCrMaxWeightNoRepeat(ia1,:);

%DataSetAllWithFirstCrMaxWeightNoRepeatCommonWeightsCathsYes = join(DataSetAllWithFirstCrMaxWeightNoRepeatCommon,exportNotNull,'Keys',{'SUBJECT_ID','ICUSTAY_ID'},'Type','inner', 'MergeKeys',true);

%  [c1,ia1,ib1] = intersect(unique(DataSetAllWithFirstCrMaxWeightNoRepeatCommon(GoodRowsindices,:).ICUSTAY_ID),unique(exportNotNull.ICUSTAY_ID1));

%  [c1,ia1,ib1] = intersect(unique(DataSetAllWithFirstCrMaxWeightNoRepeatCommon.ICUSTAY_ID),unique(exportNotNull.ICUSTAY_ID));
% find(ismember(A,B))
% [c1,ia1,ib1] = intersect(unique(DataSetAllWithFirstCrMaxWeightNoRepeat.ICUSTAY_ID),unique(urinedataNonNaN.ICUSTAY_ID));

% % [c1,ia1,ib1] = intersect(unique(DataSetAllWithFirstCrMaxWeightNoRepeatCommonWeightsCathsYes(GoodRowsindices,:).ICUSTAY_ID),unique(totalBal.ICUSTAY_ID));



%% Joining Heights
% Height2.Properties.VarNames{3} = 'HEIGHT';
% Height = union(Height1,Height2,'ICUSTAY_ID');
%Height = sortrows(Height,'ICUSTAY_ID','ascend');

%DataSetAllWithFirstCrMaxWeightNoRepeatCommonHeight = join(DataSetAllWithFirstCrMaxWeightNoRepeatCommon,Height,'Keys',{'SUBJECT_ID','ICUSTAY_ID'},'Type','leftouter', 'MergeKeys',true);

%[c1,ia1,ib1] = intersect(unique(DataSetAllWithFirstCrMaxWeightNoRepeatCommon(GoodRowsindices,:).ICUSTAY_ID),unique(Height.ICUSTAY_ID));

%% Joining Caths
% DataSetAllFirstCrMaxWeightNoReptCommonHeightWeightsCathsBoth = join(DataSetAllWithFirstCrMaxWeightNoRepeatCommonHeight,WeightsCathBoth,'Keys',{'SUBJECT_ID','ICUSTAY_ID'},'Type','inner', 'MergeKeys',true);

%%  Keep only the ones that have ICU stays in the data

%% Checking Caths
% commonCathsMaybeAndYes = CathsMaybeAndYes(ib1,:);
%
% TrueCommonCathsMaybeAndYes = length(find(commonCathsMaybeAndYes.CATH_FLAG ==1))
%
% FalseCommonCathsMaybeAndYes = length(find(commonCathsMaybeAndYes.CATH_FLAG ==0))
%
% TrueOnlyYes = length(find(DataSetAllWithFirstCrMaxWeightNoRepeatCommonWeightsCathsYes.CATH_FLAG  == 1))
%
% FalseOnlyYes = length(find(DataSetAllWithFirstCrMaxWeightNoRepeatCommonWeightsCathsYes.CATH_FLAG  == 0))

%% Getting a smaller version of the urine data.
%indicesUrineCommon = find(ismember(urinedataNonNaN.ICUSTAY_ID,DataSetAllFirstCrMaxWeightNoReptCommonHeightWeightsCathsBoth.ICUSTAY_ID));
%urinedataNonNaNCommon = urinedataNonNaN(indicesUrineCommon,:);

%% Getting a smaller version of Vasopressors

%indicesVasopressorsCommon = find(ismember(Vasopressors.ICUSTAY_ID,DataSetAllFirstCrMaxWeightNoReptCommonHeightWeightsCathsBoth.ICUSTAY_ID));
%VasopressorsCommon = Vasopressors(indicesVasopressorsCommon,:);

%% Getting a smaller version of totalBal

%indicestotalBalCommon = find(ismember(totalBal.ICUSTAY_ID,DataSetAllFirstCrMaxWeightNoReptCommonHeightWeightsCathsBoth.ICUSTAY_ID));
%totalBalCommon = totalBal(indicestotalBalCommon,:);

%% Getting a smaller version of MAP

%indicesMAPCommon = find(ismember(MAP.ICUSTAY_ID,DataSetAllFirstCrMaxWeightNoReptCommonHeightWeightsCathsBoth.ICUSTAY_ID));
%MAPCommon = MAP(indicesMAPCommon,:);

%% Change the dates to proper format.
% Vasopressors
% Dates = VasopressorsCommon.CHARTTIME;
% GoodDateFormat = FixingDateFormat(Dates);
%
% B = cell2dataset(GoodDateFormat);
% VasopressorsCommonFixed.CHARTTIME = B.CHARTTIME;
%
% % they have repeated entries for exactly the same time
%
% % MAP
% Dates = MAPCommon.CHARTTIME;
% GoodDateFormat = FixingDateFormat(Dates);

% B = cell2dataset(GoodDateFormat);
% MAPCommon.CHARTTIME = B.CHARTTIME;

%% totalBal
% totalBalCommonFixed = totalBalCommon;
% Dates = totalBalCommonFixed.CHARTTIME;
% GoodDateFormat = FixingDateFormat(Dates);
%
% B = cell2dataset(GoodDateFormat);
% totalBalCommonFixed.CHARTTIME = B.CHARTTIME;

%% Creatinine
% Dates = creatinine.CHARTTIME;
% GoodDateFormat = FixingDateFormat(Dates);
% B = cell2dataset(GoodDateFormat);
% creatinine = table2dataset(creatinine);
% creatinine.CHARTTIME = B.CHARTTIME;


%% ICUINTimes

% Dates = icuintimes.ICUSTAY_INTIME;
% GoodDateFormat = FixingDateFormat(Dates);
% B = cell2dataset(GoodDateFormat);
% icuintimes = table2dataset(icuintimes);
% icuintimes.ICUSTAY_INTIME = B.CHARTTIME;


%% Counting in in filtering UO;
global count countlength countTimes;
count = 0;% number of 'The last measurement was before the 6 hour 
%cutoff and there wasn’t any measurement after it
countlength = 0; % number of measurements < 4
countTimes = 0; 
%% Figure Settings
global width height alw fsz lw msz;
width = 3;% Width in inches
height = 3;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 18;%11;      % Fontsize
lw = 2; %1.5;      % LineWidth
msz = 8;       % MarkerSize




%% Loading data and making them global (read only never mutated so safe);


global totalBalCommonFixed MAPCommonFixed  VasopressorsCommonFixed creatinine icuintimes  urinedataNonNaNCommonTrunc;

load('matData\totalBalCommonFixed');
load('matData\MAPCommonFixed');
load('matData\VasopressorsCommonFixed');
load('matData\creatinine');
load('matData\icuintimes');
load('matData\creatinine');

load('matData\urinedataNonNaNCommon');

%% Fixing AKI creatinine definition: 0.3 or 50% increase in first 48 hours from hospital minimum.
% need to replace the current goldstandard with this.(% ToDo)

%
% [AKIcrList, TimesList,FirstCr] = ComputingAKICr2();
% save('matData\AKICrInfo.mat','AKIcrList','TimesList','FirstCr');
% %load('matData\AKICrInfo.mat');
% DataSetAllFirstCrMaxWeightNoReptCommonHeightWeightsCathsBoth.AKICrNew = AKIcrList;
% DataSetAllFirstCrMaxWeightNoReptCommonHeightWeightsCathsBoth.AKICrTruncTime = TimesList;
% DataSetAllFirstCrMaxWeightNoReptCommonHeightWeightsCathsBoth.FIRST_CR = FirstCr;
% save('matData\DataSetAllFirstCrMaxWeightNoReptCommonHeightWeightsCathsBoth.mat','DataSetAllFirstCrMaxWeightNoReptCommonHeightWeightsCathsBoth');

load('matData\DataSetAllFirstCrMaxWeightNoReptCommonHeightWeightsCathsBoth');

global data;
data = DataSetAllFirstCrMaxWeightNoReptCommonHeightWeightsCathsBoth;
clear DataSetAllFirstCrMaxWeightNoReptCommonHeightWeightsCathsBoth;


%% Truncating data
%[ urinedataNonNaNCommonTrunc ] = TruncatingUOData(); % updated 150127 with truncation based on correct creatinine trunc time
%save('matData\urinedataNonNaNCommonTrunc','urinedataNonNaNCommonTrunc');

load('matData\urinedataNonNaNCommonTrunc');
global urineDataChosen;
% 
% 'number of patients removed from truncated version of urine output'
% length(unique(urinedataNonNaNCommon.ICUSTAY_ID)) - length(unique(urinedataNonNaNCommonTrunc.ICUSTAY_ID))
whichUrine = 'Trunc'; % 'Reg' or 'Trunc'

if strcmp(whichUrine, 'Trunc')
    
    urineDataChosen = urinedataNonNaNCommonTrunc;
    clear urinedataNonNaNCommon urinedataNonNaNCommonTrunc;
    
elseif strcmp(whichUrine,'Reg')
    urineDataChosen = urinedataNonNaNCommon;
    clear urinedataNonNaNCommon urinedataNonNaNCommonTrunc;
    
end


%% PreAnalyses (produces histograms): not used much: to 

%[ numMeasAll,numMeas6hoursAll, timeIntervalsAll, timeIntervals6hoursAll] = filterAndPreprocessing_analysis();
%save('matData\AfterfilterAndPreprocessing_analysis.mat','numMeasAll','numMeas6hoursAll', 'timeIntervalsAll', 'timeIntervals6hoursAll');
load('matData\AfterfilterAndPreprocessing_analysis')


%% Processing and filtering
%[GoodRowsindices,GoodRowsTF,Values]  = filterAndPreprocessing();% GoodRowsTF is to be used with other ones NaN
%save('matData\AfterfilterAndPreprocessing.mat','GoodRowsindices','GoodRowsTF','Values');
load('matData\AfterfilterAndPreprocessing.mat') 

Histograms(Values);
%% estimating LBMs: one time only because saving to dataset.

%  [ EstimatedLBMs ] = EstimateLBM();
%  DataSetAllFirstCrMaxWeightNoReptCommonHeightWeightsCathsBoth.LBM = EstimatedLBMs;
%  save('DataSetAllFirstCrMaxWeightNoReptCommonHeightWeightsCathsBoth.mat','DataSetAllFirstCrMaxWeightNoReptCommonHeightWeightsCathsBoth');
% Note: not having filter on control group

%% Calculate goldstandard and make binary and check NaN's

%% choose between override or not blocks%%

%% with override
% AKI_OVERIDEind = 13;
% AKI_CATEGORY_ICU48ind = 5;% using hospital min
%
% AKI_OVERIDE_COl = ConvertToCell(data,AKI_OVERIDEind,AllRows);
%
% AKI_CATEGORY_ICU48_COl = ConvertToCell(data,AKI_CATEGORY_ICU48ind,AllRows); % hospital min
%
% GoldStandardColumn = computeGoldStandard(AKI_OVERIDE_COl,AKI_CATEGORY_ICU48_COl);% is of cell array type
% GoldStandardColumnMat = cell2mat(GoldStandardColumn);




%% withoutoverride uncomment if used
AllRows = 1: size(data,1);
AKI_CATEGORY_ICU48ind = 75;% using hospital min  % new AKI

GoldStandardColumn = ConvertToCell(data,AKI_CATEGORY_ICU48ind,AllRows); % hospital min

GoldStandardColumnMat = cell2mat(GoldStandardColumn);

% see which of the goldstandard rows ~= NaN

TFNonNullGoldstandards =  ~isnan(GoldStandardColumnMat);

% converging categories 1 ,2  and 3 into 1; Not needed anymore
NonNullGoldStandardRows = find(TFNonNullGoldstandards==1);
%GoldStandardColumnMat(NonNullGoldStandardRows,1) = GoldStandardColumnMat(NonNullGoldStandardRows,1)>=1;

% update good indices
GoodRowsTF = GoodRowsTF & TFNonNullGoldstandards;
GoodRowsindices = find(GoodRowsTF ==1);

'No Gold Standard'
length(AllRows) - length(NonNullGoldStandardRows);
%% Checking if change in definition made many differences; using previous def as gold gold
% 
% AKIcrList
% GoldStandardColumnMat
% [CM,order] = confusionmat(AKIcrList,GoldStandardColumnMat,'ORDER',[1,0])
% ConfusionMatrixInfo(CM);
% AKIcrList
% TFAKIcrList =  ~isnan(AKIcrList);
% 
% NonNullGoldStandardRowsNew = find(TFAKIcrList==1);

% %% Bad rows ICU stay histogram : don't run this block
% BadRows = setdiff(AllRows,GoodRowsindices);
% ICUColumn = ConvertToCell(data,14,BadRows);
% 
% % CSRU
% TFMatchesICU = (strcmp('CSRU',ICUColumn));
% RowsMatchingICU = find(TFMatchesICU==1);
% numCSRU = length(RowsMatchingICU);
% 
% % SICU
% TFMatchesICU = (strcmp('SICU',ICUColumn));
% RowsMatchingICU = find(TFMatchesICU==1);
% numSICU = length(RowsMatchingICU);
% 
% % CCU
% TFMatchesICU = (strcmp('CCU',ICUColumn));
% RowsMatchingICU = find(TFMatchesICU==1);
% numCCU = length(RowsMatchingICU);
% 
% %MICU
% ICU1 = 'MICU';% 'MICU'+'FICU'
% ICU2 = 'FICU';
% TFMatchesICU = (strcmp(ICU1,ICUColumn) | strcmp(ICU2,ICUColumn));
% RowsMatchingICU = find(TFMatchesICU==1);
% numMICU = length(RowsMatchingICU);
% 
% 
% figure()
% pos = get(gcf, 'Position');
% set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
% set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
% bar([numCSRU,numSICU,numCCU,numMICU]);
% set(gca,'xTickLabel',{'CSRU','SICU','CCU','MICU'});
% title('ICUSTAY FIRST SERVICE - excluded');
% print('ICUSTAY FIRST SERVICE - excluded','-dpng','-r300');
% 
% close all


%% Checking ethnicity.

[idx,label] = grp2idx(ConvertToCell(data,16,GoodRowsindices)); % AllRows
NumberPerEth = NaN*ones(40,1);

for ethIdx = 1: 40 % might change this
    NumberPerEth(ethIdx) = length(find(idx == ethIdx));
end


%% Table 1 analyses using GoodRowsindices based on only the filtering - note this is not the table 1 that is used in the manuscript, 
% there is another table that is generated using CreateTable1_v3 that is
% after the GEE part as it needs the GEE cohort.

% which ICU  % All, CSRU, 	SICU,	CCU,	MICU + FICU

%Table1 = CreateTable1( data, AllRows, GoodRowsTF);
%save('matData\Table1.mat','Table1');
load('matData\Table1.mat');

% Print count countlength countTimes;
'count',count
'countlength',countlength
'countTimes',countTimes
%% GetFeatures


ICUSTAY_ID_COL = cell2mat(ConvertToCell(data, 2,AllRows));

AllrowsUrine = 1: size(urineDataChosen,1);
UrineICU = cell2mat(ConvertToCell(urineDataChosen, 2,AllrowsUrine)); % from urine data, different than in data
TimeThresholdRange = [2,4,6,8,10,12,14,16,18,20,22,24];%[1,2,3,4,5,6,7,8,9,10,11,12];%[6];%[5,6,7];%[3,4,5,6,7,8,9];%[3:1:24];%[3:1:12];   % [3:24] : colSub
global AvgVolumeThresholdRange;
AvgVolumeThresholdRange = [0: 0.1: 1];
numIndices = length(TimeThresholdRange) * length(AvgVolumeThresholdRange);


%%         %%%%%%%%%%% Only dealing with GoodRowsindices %%%%%%%%%%%%%%%%%%%%
ICUSTAY_ID_COL = ICUSTAY_ID_COL(GoodRowsindices,1);
GoldStandardColumnMat = GoldStandardColumnMat(GoodRowsindices,1);

AKIUOCurrentdef  =  NaN*ones(size(GoodRowsindices,1),numIndices);
AKIUOCurrentdefVol = AKIUOCurrentdef;

NewData = [ICUSTAY_ID_COL GoldStandardColumnMat AKIUOCurrentdef]; % true false, AKI per comb of (Tth,Vth)
NewDataVols = [ICUSTAY_ID_COL GoldStandardColumnMat AKIUOCurrentdefVol]; % min volumes, would be repeated for Vth to be consistent
NewDataStartTimeMins = [ICUSTAY_ID_COL GoldStandardColumnMat AKIUOCurrentdefVol];
%

% controls
totalFluidPriorMin = [ICUSTAY_ID_COL GoldStandardColumnMat AKIUOCurrentdefVol];

TFVasopressorPriorBelowMin = [ICUSTAY_ID_COL GoldStandardColumnMat AKIUOCurrentdefVol];

avgMAP1hrPriorBelowMin = [ICUSTAY_ID_COL GoldStandardColumnMat AKIUOCurrentdefVol];

% for each of the time window lengths there is a column
[ ControlVars] = ConvertControlVarsToMat(AllRows,data);
ControlVarsGood = ControlVars(GoodRowsindices,2:end) ; % not using ICU stay now
% ControlVars is based on original indices ;
% ControlVarsGood is based on GoodIndices indicization
ForPreAllocation = NaN*ones(100* 11050,1);
ForPreAllocationX = NaN*ones(100* 11050,12);

GoldGEE = cell(1,length(TimeThresholdRange));  
IDGEE = cell(1,length(TimeThresholdRange)); % use ICUStayID
TimesGEE = cell(1,length(TimeThresholdRange)); % use in seconds from first measurement
XGEE = cell(1,length(TimeThresholdRange)); % Control Vars, with Volume
% 'AvgVolumePerwindow','Gender','Age','Weight', 'Diuretics', 'Height', 'FirstCr', 'Diabetes',
% 'LBM', 'avgMAP1hrPriorWindow', 'totalFluidPriorWindow', 'TFVasopressorsPriorWindow'

% preallocating the nested cell arrays
for timeIndx = 1:length(TimeThresholdRange)
    
    GoldGEE{timeIndx}= ForPreAllocation;
    IDGEE{timeIndx} = ForPreAllocation;
    TimesGEE{timeIndx} = ForPreAllocation;
    XGEE{timeIndx} = ForPreAllocationX;
end

IndexForArrays = ones(1,length(TimeThresholdRange));

% Making Function handles to be able to later vectorize computations and
% make it more efficient
ComputeTotalFluidPriorFuncHandle = @(IcuID,StartTime, firstTimeStr) ComputeTotalFluidPrior(IcuID,StartTime, firstTimeStr);
ComputeTFVasopressorPriorFuncHandle =  @(IcuID,StartTime, firstTimeStr)ComputeTFVasopressorPrior(IcuID,StartTime, firstTimeStr);
ComputeAvgMAP1hrPriorFuncHandle =  @(IcuID,StartTime, firstTimeStr) ComputeAvgMAP1hrPrior(IcuID,StartTime, firstTimeStr);
ComputeAKIGiven4FuncHandle = @(Subdata,TimeThreshold,FirstWeightValue)ComputeAKIGiven4 (Subdata,TimeThreshold,FirstWeightValue);         


close all;
n = 1; % change back 
save('matData\beforeBigLoop.mat');
%load('beforeBigLoop.mat');
for IcuID = ICUSTAY_ID_COL' %unique(ICUSTAY_ID_COL')
    % find all the rows of urineDataChosen that have the ICUID
    
    FoundICUIdxs = find(UrineICU(:,1) == IcuID ) ;
    
    if ~isempty(FoundICUIdxs)
        % Subdata
        Subdata = urineDataChosen(FoundICUIdxs,:);
        
        OrigIndex = GoodRowsindices(n);
        FirstWeightValue = data(OrigIndex,68).FIRST_WEIGHT;
        
        
        % returns a list for AKIUOCurrentdefPerCombfinal and
        % StartTimeBelowThreshold instead of just a number
        
        
        repSubdataCell = repcell(Subdata,length(TimeThresholdRange),1);% should be column
        
        TimeThresholdRangeCell = mat2cell(TimeThresholdRange',ones(1,length(TimeThresholdRange),1));
        
        repFirstWeightValueCell = repcell(FirstWeightValue,length(TimeThresholdRange),1);
        
        % for each patient see if they have AKI based on the different UO criteria
        [AKIUOCurrentdefPerCombfinalList_AllTimes,minVol_AllTimes, startTimeMin_AllTimes, firstTimeStr_AllTimes,lastTimeStr_AllTimes, AKIUOCurrentdefPerCombListVol_AllTimes,StartTimes_AllTimes] =  cellfun(ComputeAKIGiven4FuncHandle,repSubdataCell,TimeThresholdRangeCell,repFirstWeightValueCell,'UniformOutput',false);% stopped here
        
        
        repICUminAllTimescell = repcell(IcuID,length(TimeThresholdRange),1);
        
        
        avgMAP1hrPriorBelowMinVal_AllTimes = cellfun(ComputeAvgMAP1hrPriorFuncHandle,repICUminAllTimescell,startTimeMin_AllTimes, firstTimeStr_AllTimes,'UniformOutput',false);
        TFVasopressorPriorBelowMinVal_AllTimes = cellfun(ComputeTFVasopressorPriorFuncHandle,repICUminAllTimescell,startTimeMin_AllTimes, firstTimeStr_AllTimes,'UniformOutput',false);
        totalFluidPriorMinVal_AllTimes = cellfun(ComputeTotalFluidPriorFuncHandle,repICUminAllTimescell,startTimeMin_AllTimes, firstTimeStr_AllTimes,'UniformOutput',false);
        
        
        % colSub = indx for time threshold. rowSub = indx for volume threshold
        %rowSub =0;
        colSub =0;
        for TimeThreshold = TimeThresholdRange
            colSub =colSub +1;
            
            linearInds = sub2ind([length(AvgVolumeThresholdRange),length(TimeThresholdRange)], 1:length(AvgVolumeThresholdRange), repmat(colSub, 1, length(AvgVolumeThresholdRange)));
            
            NewData(n,linearInds+2)= AKIUOCurrentdefPerCombfinalList_AllTimes{colSub};%AKIUOCurrentdefPerCombfinalList;% now doing all volumeThresholds
            
            NewDataVols(n,linearInds+2) = minVol_AllTimes{colSub};% used in decision tree.
            
            NewDataStartTimeMins(n, linearInds+2) = startTimeMin_AllTimes{colSub};
            
            totalFluidPriorMin(n,linearInds+2) = totalFluidPriorMinVal_AllTimes{colSub};
            TFVasopressorPriorBelowMin(n,linearInds+2) = TFVasopressorPriorBelowMinVal_AllTimes{colSub};
            avgMAP1hrPriorBelowMin(n,linearInds+2) = avgMAP1hrPriorBelowMinVal_AllTimes{colSub};
            
            
            %% Going through each window and extracting for GEE/
            % replicating values
            
            
            GoldGEE{colSub}(IndexForArrays(colSub):IndexForArrays(colSub)+length(StartTimes_AllTimes{colSub})-1)  = repmat(GoldStandardColumnMat(n), length(StartTimes_AllTimes{colSub}), 1);
            IDGEE{colSub}(IndexForArrays(colSub):IndexForArrays(colSub)+length(StartTimes_AllTimes{colSub})-1) = repmat(IcuID, length(StartTimes_AllTimes{colSub}), 1) ;
            TimesGEE{colSub}(IndexForArrays(colSub):IndexForArrays(colSub)+length(StartTimes_AllTimes{colSub})-1) = StartTimes_AllTimes{colSub}' ;
            
            
            
            repICUcell = repcell(IcuID,length(StartTimes_AllTimes{colSub}),1);
            repfirstTimeStrCell = repcell(firstTimeStr_AllTimes{colSub},length(StartTimes_AllTimes{colSub}),1);
            StartTimesCell = mat2cell(StartTimes_AllTimes{colSub}',ones(1,length(StartTimes_AllTimes{colSub}),1));
            
            TotalFluidPriorWindows = cellfun(ComputeTotalFluidPriorFuncHandle,repICUcell, StartTimesCell,repfirstTimeStrCell);
            TFVasopressorPriorWindows =  cellfun(ComputeTFVasopressorPriorFuncHandle,repICUcell, StartTimesCell,repfirstTimeStrCell);
            AvgMAP1hrPriorWindows =  cellfun(ComputeAvgMAP1hrPriorFuncHandle,repICUcell, StartTimesCell,repfirstTimeStrCell);
            
            X = [AKIUOCurrentdefPerCombListVol_AllTimes{colSub}' repmat(ControlVarsGood(n,:),length(StartTimes_AllTimes{colSub}),1) TotalFluidPriorWindows TFVasopressorPriorWindows  AvgMAP1hrPriorWindows]; % Control Vars
            
            XGEE{colSub}(IndexForArrays(colSub):IndexForArrays(colSub)+length(StartTimes_AllTimes{colSub})-1,:) = X;
            % 'avgVol', 'Gender','Age','Weight', 'Diuretics', 'Height', 'FirstCr', 'Diabetes',
            % 'LBM', 'totalFluidPriorWindow','TFVasopressorsPriorWindow''avgMAP1hrPriorWindow', 
            
            IndexForArrays(colSub) = IndexForArrays(colSub) + length(StartTimes_AllTimes{colSub});
        end
        
    end
    
    if mod(n,1000) ==0
        n
        save('matData\DuringBigLoop.mat','NewData','NewDataVols','NewDataStartTimeMins', 'totalFluidPriorMin',  'TFVasopressorPriorBelowMin',  'avgMAP1hrPriorBelowMin','GoldGEE','IDGEE','TimesGEE','XGEE','IndexForArrays');

    end
    if n < 5
        n
    end
    n = n+1;
end
 
%save('matData\AfterBigLoop.mat','NewData','NewDataVols','NewDataStartTimeMins', 'totalFluidPriorMin',  'TFVasopressorPriorBelowMin',  'avgMAP1hrPriorBelowMin','GoldGEE','IDGEE','TimesGEE','XGEE','IndexForArrays');
%save('matData\EverythingAfterBigLoop.mat');
load('matData\EverythingAfterBigLoop.mat')
% Note: the values are in seconds

%load('matData\EverythingAfterBigLoop.mat');
%% Converting inf to NaN; because INFs causes the matrix to
% become singular.
[I,J] = find(NewDataVols(:,:) == inf); NewDataVols(I,J)= NaN*ones(size(I,1),size(J,1)); 

%% GEE preprocessing

% label Cr meas and truncate
TruncateLabelCreatinineMeas() % saves creatinineTruncLabel
IDGEEOrig = IDGEE; GoldGEEOrig = GoldGEE; TimesGEEOrig = TimesGEE; XGEEOrig =XGEE;
% Remove crossing windows, make time GEE relative to admission, 
[changes,IDGEE,GoldGEE,TimesGEE,XGEE ] = FilteringCrossWindAndLabeling(IDGEEOrig,GoldGEEOrig, TimesGEEOrig, XGEEOrig );
%save('GEEVars_ noCrossingWindows.mat','IDGEE','GoldGEE','TimesGEE','XGEE');
load('GEEVars_ noCrossingWindows.mat');
% If needed;
%  w = warning('query','last');
% id = w.identifier;
% warning('off',id)
%% Adding (In - Out)/Weight;
%10:'totalFluidPriorWindow'; 1:'avgVol'

XGEEbefore = XGEE;

for timeDur = 1:12
    timeDurationValue = TimeThresholdRange(timeDur);
    currentX = XGEE{timeDur};
    % out[vol] = out[vol/kg*hrs] *weight*time duration
    Out_vol = currentX(:,1).*currentX(:,4)*timeDurationValue;
    
    In_Out_Weight = (currentX(:,10) - Out_vol)./currentX(:,4);
    XGEE{timeDur} = [currentX In_Out_Weight];
    
end


%% TF CCU
ICUID_allrows = data.ICUSTAY_ID;
ICUGEE = cell(1,12); % 1 : CSRU,  2: SICU, 3: CCU,  4:MICU + FICU

for timeDur = 1:12
    ICUGEE{timeDur} = NaN*ones(length(IDGEE{timeDur}),1);
    IDGEE_time = IDGEE{timeDur};
    [lia1,locb1] = ismember(IDGEE_time,ICUID_allrows); % None of them should not be found.
    ICUGEE{timeDur} = ControlVars(locb1,1);
       
end

% %% Adding fluid out/fluid in to XGEE
% %10:'totalFluidPriorWindow'; 1:'avgVol'
% XGEEbefore = XGEE;
% 
% for timeDur = 1:12
%     currentX = XGEE{timeDur};
%     XGEE{timeDur} = [currentX currentX(:,1)./currentX(:,10)];
%     
% end

%% GEE analyses
addpath('GEEQBOX'); % from: get website link : ToDo

% warning using the same variable names as ROC analyses this might cause
% some error if not careful: when cleaning fix this. 

%load('GEEVars_ noCrossingWindows.mat');

% % 'avgVol', 'Gender','Age','Weight', 'Diuretics', 'Height', 'FirstCr', 'Diabetes', % 'LBM', 'totalFluidPriorWindow','TFVasopressorsPriorWindow''avgMAP1hrPriorWindow', 
%[betahat,alphahat,results] = gee(id,percent,month,X,'b','markov',varnames);
varNamesOrig = {'avgVol', 'Gender','Age','Weight', 'Diuretics', 'Height', 'FirstCr', 'Diabetes','LBM','totalFluidPriorWindow','TFVasopressorsPriorWindow','avgMAP1hrPriorWindow','In_Out_Weight','cons'};
%1:'avgVol', 2:'Gender',3: 'Age',4:'Weight', 5:'Diuretics', 6:'Height', 7:'FirstCr', 8:'Diabetes',9:'LBM',10:'totalFluidPriorWindow',11:'TFVasopressorsPriorWindow',12:'avgMAP1hrPriorWindow',13:'In_Out_Weight',14:'cons'
ToKeepVarsInd = [3,5,7,9,11:13]; %[2:9,11:13];%[1:13];%[1,3,5,7,8,9];
ToKeep = [ToKeepVarsInd,14]; % cons should always be there but XGEE doesn't have cons column
varNames = varNamesOrig(ToKeep);


%[betahat,alphahat,results] = gee(id,percent,month,X,'b','markov',varnames);
% cons should be the last variable
%varNames = {'avgVol','cons'};%, 'Gender','Age','Weight', 'Diuretics', 'Height', 'FirstCr', 'Diabetes','LBM','totalFluidPriorWindow','TFVasopressorsPriorWindow','avgMAP1hrPriorWindow','cons'};


models = cell(1+length(TimeThresholdRange),1+ length(varNames)+2);
models{1,1} = 'Time dur';
models(1,2:end-2) = varNames;
models(1,end-1:end)= {'AUC','distJPoint'};
%'GoldGEE','IDGEE','TimesGEE','XGEE'



UniqueIDs = cell(1,12); 
UniqueIDsUsed = cell(1,12); 
numWindows = NaN*ones(1,12);
AUCvals = NaN*ones(1,12);
CIhighMat = NaN*ones(12,length(ToKeep));
CIlowMat = NaN*ones(12,length(ToKeep));
CoeffEstMat = NaN*ones(12,length(ToKeep));
pValsMat = NaN*ones(12,length(ToKeep));

ChooseICU = 0; % 0: do not segregate, 1: do segregate
ICUchoiceTF = 0; % choose ICU == 3: CCU


% for comparing to 6, 0.5 thrshold prediction 

Ind_pred__6_0_5 =  2 + sub2ind([length(AvgVolumeThresholdRange),length(TimeThresholdRange)], find(AvgVolumeThresholdRange == 0.5), find(TimeThresholdRange == 6));
ICUID_Gold_Pred__6_0_5_Table =  array2table(NewData(:,[1,2,Ind_pred__6_0_5]), 'VariableNames', {'ICUSTAY_ID', 'AKIGold', 'Pred__6_0_5'});
NRI_IDI_GEE_models_Table_colnames = {'TimeThreshold', 'NRI', 'NRI_pval', 'IDI', 'IDI_pval'};
NRI_IDI_GEE_models_Table = array2table(nan*ones(length(TimeThresholdRange),length(NRI_IDI_GEE_models_Table_colnames)), 'VariableNames', NRI_IDI_GEE_models_Table_colnames);




randomSeedNumList = 1; %1:10; % used 1,2,3  % if just want to 

for randomSeedNum = randomSeedNumList
    rng(randomSeedNum) % set seed for randomizer  in selecting the training and test sets in GEE models 

    mkdir('GEE_NRI_IDI')
    for colSub = 1: length(TimeThresholdRange) % 1: length(TimeThresholdRange)

        models{1+colSub,1} = TimeThresholdRange(colSub); 

        % extracting GEE variables for the specified time duration

        GoldGEETimeNotFiltered = GoldGEE{colSub} ;
        IDGEETimeNotFiltered = IDGEE{colSub};
        TimesGEETimeNotFiltered = TimesGEE{colSub};
        XGEETimeNotFilteredAllFeat = XGEE{colSub};
        ICUGEETimeNotFiltered = ICUGEE{colSub};


        XGEETimeNotFiltered = XGEETimeNotFilteredAllFeat(:,ToKeepVarsInd);%%% Added 150209 to not remove rows based on vars not used 



    %     % Remove NaN because not filled till that point
    %     [IsNaNID,JsNaNID] = find(isnan(IDGEETime)==1);
    %     GoldGEETime(IsNaNID) = [];
    %     IDGEETime(IsNaNID) = [];
    %     TimesGEETime(IsNaNID) = [];
    %    
    %     XGEETime(IsNaNID,:) = [];
    %     ICUGEETime(IsNaNID) = [];


        % Remove NaN and all associated rows in all GEE variables

        [IsNaN,JsNaN] = find(isnan(XGEETimeNotFiltered)==1);
         GoldGEETime  = GoldGEETimeNotFiltered; IDGEETime = IDGEETimeNotFiltered; TimesGEETime = TimesGEETimeNotFiltered; XGEETime = XGEETimeNotFiltered;ICUGEETime = ICUGEETimeNotFiltered;
        
        GoldGEETime(IsNaN) = [];
        IDGEETime(IsNaN) = [];
        TimesGEETime(IsNaN) = [];

        XGEETime(IsNaN,:) = [];
        ICUGEETime(IsNaN,:) = [];

        % Select only people in CCU or not? :
        if ChooseICU
            if ICUchoiceTF  == 1 % select only ones in CCU
                rowsInICU = find(ICUGEETime == 3);
            end
            if ICUchoiceTF == 0 
                rowsInICU = find(ICUGEETime ~= 3);    
            end

            GoldGEETime = GoldGEETime(rowsInICU);
            IDGEETime = IDGEETime(rowsInICU);
            TimesGEETime = TimesGEETime(rowsInICU);

            XGEETime = XGEETime(rowsInICU,:);
            ICUGEETime = ICUGEETime(rowsInICU);


        end


        % if var 13:  included: In_Out_Weight  (fluid balance = (fluid
        % input - fluid output )/ weight)  % sometimes weight ~ 0 which
        % makes fluid balance = Inf 
        
         [LIA,LOCB] = ismember(13,ToKeep);
         if LIA
             [Iinf,Jinf] = find(isinf(XGEETime(:,LOCB)));
             GoldGEETime(Iinf) = [];
             IDGEETime(Iinf) = [];
             TimesGEETime(Iinf) = [];

             XGEETime(Iinf,:) = [];

        end
        % saving uniqueIDs and number of windows for each time duration
        UniqueIDs{colSub} = unique(IDGEETime);
        numWindows(colSub) = length(IDGEETime);


        % Getting Training and Validation indices
        % doing this because qls complains that it is out of memory.
        maxNumber = 80000;%167192;
        %maxNumber = ceil(55000 + 0.3*(10/7)*55000);

        if length(GoldGEETime) < maxNumber % tested various amounts in qls.
            [trainInd,valInd,testInd] = dividerand(length(GoldGEETime),.7,0.3,0);

        else
            strcat('Using maxNumber, for time duration:  ' ,  num2str(TimeThresholdRange(colSub)))
            [trainInd,valInd,testInd] = dividerand(maxNumber,.7,0.3,0);
            UniqueIDsUsed{colSub} = unique(IDGEETime(1:maxNumber));
        end
 
        
        
        
        GoldGEETimeTraining = GoldGEETime(trainInd);
        GoldGEETimeValid =  GoldGEETime(valInd);

        IDGEETimeTraining = IDGEETime(trainInd);
        IDGEETimeValid = IDGEETime(valInd);

        TimesGEETimeTraining = TimesGEETime(trainInd);
        TimesGEETimeIDGEEValid = TimesGEETime(valInd);

        XGEETimeTraining = [XGEETime(trainInd,:) ones(length(trainInd),1)];% changed 
        XGEETimeValid = [XGEETime(valInd,:) ones(length(valInd),1)]; % changed

        % getting model
        % need Beta values and 95% CI values for each.
        %tillInd = 60000; % works: 50000, 70000,75000, 78800
        %[betahat,alphahat,results] = qls(IDGEETimeTraining(1:tillInd),GoldGEETimeTraining(1:tillInd),TimesGEETimeTraining(1:tillInd),XGEETimeTraining(1:tillInd,:),'b','markov',varNames); 
        [betahat,alphahat,results] = qls(IDGEETimeTraining,GoldGEETimeTraining,TimesGEETimeTraining,XGEETimeTraining,'b','markov',varNames); % changed here

        %(1:33260)
        CILowUp =  results.robust(3:end,6:7); % from robust
        pValue = results.robust(3:end,5);

        %PredictedAKI = sum(bsxfun(@times,XGEETimeValid,betahat'),2); % column multiply
        %1 / ( 1 + exp(-X) )
        PredictedAKI = 1./(1+exp(-1*sum(bsxfun(@times,XGEETimeValid,betahat'),2))); % column multiply; rows are not each ICU ID, they are each window


        [X,Y,T,AUC, OPTROCPT] = perfcurve(GoldGEETimeValid,PredictedAKI,[1]); %% Probabilty or True, T = thresholds, X, Y : 1- spec , sens?
        [distJpoint,minInd] = nanmin(sqrt((X-0).^2  + (Y-1).^2)); % note: the Jpoint 


        %% Perform NRI, IDIs tatistics comparing to the 6, 0.5 threshold
        % definition
        JpointProbThresh = T(minInd); % find probability threshold for J point
        binPredAKIValidTime = PredictedAKI >= JpointProbThresh; %convert the predictions of the test data to binary - still at window level



        % In order to compare to 6, 0.5 need to extract predictions of AKI from window level to patient level
        ICUID_NewPredTable_windLevel = array2table([IDGEETimeValid binPredAKIValidTime],'VariableNames', {'ICUSTAY_ID','binPredAKIValidTime'});
        ICUID_NewPred_ICUIDLevel = grpstats(ICUID_NewPredTable_windLevel, 'ICUSTAY_ID', 'nanmax'); % like group_by and summarize in R

        ICUID_Gold_Pred__6_0_5__NewPredTable = innerjoin(ICUID_NewPred_ICUIDLevel,ICUID_Gold_Pred__6_0_5_Table,'Keys',{'ICUSTAY_ID'});

        % write table of predictions -  not a result  - just used to see if R
        % method would give the same answer
       % writetable(ICUID_Gold_Pred__6_0_5__NewPredTable, "GEE_NRI_IDI\ICUID_Gold_Pred__6_0_5__NewPredTable" +num2str(TimeThresholdRange(colSub)) + ".txt",'WriteRowNames',true, 'Delimiter','\t')  


        % comparing to original threshold of 6, 0.5
        % NRI

        NRI_IDI_GEE_models_Table.TimeThreshold(colSub) = TimeThresholdRange(colSub);
        [ NRI , NRI_pval ] = NetReclassificationImprovement(ICUID_Gold_Pred__6_0_5__NewPredTable.Pred__6_0_5, ICUID_Gold_Pred__6_0_5__NewPredTable.nanmax_binPredAKIValidTime, ICUID_Gold_Pred__6_0_5__NewPredTable.AKIGold); % [ NRI , pval ] = NetReclassificationImprovement( pred_old , pred_new , outcome )
        [ IDI , IDI_pval ] = IntegratedDiscriminationImprovement(ICUID_Gold_Pred__6_0_5__NewPredTable.Pred__6_0_5, ICUID_Gold_Pred__6_0_5__NewPredTable.nanmax_binPredAKIValidTime, ICUID_Gold_Pred__6_0_5__NewPredTable.AKIGold); % [ IDI , pval ] = IntegratedDiscriminationImprovement( pred_old , pred_new , outcome )

         NRI_IDI_GEE_models_Table.NRI(colSub) = NRI;
         NRI_IDI_GEE_models_Table.NRI_pval(colSub) = NRI_pval;

         NRI_IDI_GEE_models_Table.IDI(colSub) = IDI;
         NRI_IDI_GEE_models_Table.IDI_pval(colSub) = IDI_pval;


        %% extract coefficients and other statistics for each of the features
        for controlInd = 1:length(varNames)
            models{1+colSub,1+controlInd} = ['(',num2str(CILowUp{controlInd,1}),', ',num2str(betahat(controlInd)),', ', num2str(CILowUp{controlInd,2}),', ',num2str(pValue{controlInd}), ')'];

            CIlowMat(colSub,controlInd) = CILowUp{controlInd,1};
            CoeffEstMat(colSub,controlInd) = betahat(controlInd);
            CIhighMat(colSub,controlInd) = CILowUp{controlInd,2};
            pValsMat(colSub,controlInd) = pValue{controlInd};

        end

        models{1+colSub, end-1} = AUC;
        models{1+colSub, end} = ['(',num2str(distJpoint),', ',num2str(X(minInd)),', ',num2str(Y(minInd))];

        AUCvals(colSub) = AUC;

        figure();
        plot(X,Y, 'o');
        xlabel('False positive rate'); ylabel('True positive rate')
        figName = ['ROC - GEE time duration ',num2str(TimeThresholdRange(colSub))];
        title(figName);
        text([X(1) X(end)],[Y(1) Y(end)],{num2str(T(1)),num2str(T(end))}) % label the first last and best cutoff
        text(X(minInd),Y(minInd),num2str(T(minInd)));
        print(figName,'-dpng','-r300');

        close;
    end
   % writetable(NRI_IDI_GEE_models_Table, "GEE_NRI_IDI\NRI_IDI_GEE_models_Table_All_FludBal_noInOut_NoWHG_rng_" + num2str(randomSeedNum)  + ".xlsx", 'WriteRowNames',true)  


end

modelsTable = cell2table(models);
writetable(modelsTable , "GEEResults\modelsTable.xlsx"); % saves as header model1, model2, ..; need to open the excel file and remove the first row


% convert to table, add column time threshold
CIhighTable = array2table(CIhighMat, 'VariableNames', varNames); CIhighTable.TimeThreshold = TimeThresholdRange';
CIlowTable = array2table(CIlowMat, 'VariableNames', varNames); CIlowTable.TimeThreshold = TimeThresholdRange';
CoeffEstTable = array2table(CoeffEstMat, 'VariableNames', varNames); CoeffEstTable.TimeThreshold = TimeThresholdRange';
pValsTable = array2table(pValsMat,  'VariableNames', varNames); pValsTable.TimeThreshold = TimeThresholdRange';

AUCvalsTable = table(); AUCvalsTable.AUCVals = AUCvals'; AUCvalsTable.TimeThreshold = TimeThresholdRange';

mkdir('GEEResults');
writetable(CIhighTable, "GEEResults\CIhighTable.xlsx"); writetable(CIlowTable, "GEEResults\CIlowTable.xlsx"); 
writetable(CoeffEstTable, "GEEResults\CoeffEstTable.xlsx"); writetable(pValsTable, "GEEResults\pValsTable.xlsx"); 
writetable(AUCvalsTable, "GEEResults\AUCvalsTable.xlsx"); 

save('matData\variableCoeffsForGEEModels__All_FludBal_noInOut_NoWHG.mat', 'CIhighTable', 'CIlowTable', 'CoeffEstTable', 'pValsTable', 'AUCvalsTable');

%save('matData\models_All_FludBal_noInOut_NoWHG.mat','models','UniqueIDsUsed','UniqueIDs','numWindows','AUCvals');
load('matData\models_All_FludBal_noInOut_NoWHG.mat');  % is the one used in paper
%save('matData\modelsAvgVolOnly.mat','models','UniqueIDsUsed','UniqueIDs');
%save('matData\modelsAll.mat','models')

%% Analyse the NRI for the GEE models for various random number seeds
analyseNRI_GEE(); 




%% Generating Table 1

%save('matData\models.mat','models','UniqueIDsUsed','UniqueIDs');
cellsz = cellfun(@size,UniqueIDsUsed,'uni',false);
a =  UniqueIDs;
dim =1 ;
result = sum(cellfun(@(c) size(c,dim), a), dim);
%% Table 1 modified with the ones entering GEE
chosenTimeInd = 8; % Tth = 16

ICUIDsusedGEE = union(UniqueIDsUsed{chosenTimeInd},UniqueIDs{1,chosenTimeInd} );
GEEgoodrowTF = ismember(data.ICUSTAY_ID,ICUIDsusedGEE);
%Table1_mod = CreateTable1_v2(GoodRowsTF, GEEgoodrowTF);
Table1_mod_3 = CreateTable1_v3(GoodRowsTF, GEEgoodrowTF);




%% Histogram of first Cr for the people that entered GEE

varIndx = 7;
%7:'FirstCr'
for colSub = 1:12
    fig = figure();
    
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
    set(gca, 'FontSize', fsz);
    
    timeDur = TimeThresholdRange(colSub);
    XGEETime = XGEE{colSub};
    FirstCrTime = XGEETime(:,varIndx); 
    minFirstCr = nanmin(FirstCrTime);
    maxFirstCr = nanmax(FirstCrTime);
    medianFirstCr = nanmedian(FirstCrTime);
        
    hist(FirstCrTime,50);
    
    
    xlabel('First creatinine value')
    ylabel('Count')
    figName = ['Histogram of first creatinine value - timeDur',num2str(timeDur)];
    title([figName,' - min ',num2str(minFirstCr) ,' max ', num2str(maxFirstCr) ,' median ',num2str(medianFirstCr) ]);
    print(figName,'-dpng','-r300');
    close;
end




%% ROC Curves summaries  for each window size and for each control variable - using mnrfit - multinomial regression for nominal responses - not used in paper
% Splitting data into training and validation sets
% rows = different time durations
% columns = different control variables

%save('TillNowBeforeROC_FixedAKIUONotTrunc.mat');
save('TillNowBeforeROC_Trunc.mat');

load('TillNowBeforeROC_FixedAKIUONotTrunc.mat');
%load('TillNowBeforeROC_Trunc.mat');
AUCTable = cell(length(TimeThresholdRange), 13 );

            
ControlVarNames = {'Time Duration', 'No ctrl','Gender','Age','Weight', 'Diuretics', 'Height', 'FirstCr', 'Diabetes','LBM', 'totalFluidPriorMin', 'TFVasopressorsPriorMin','avgMAP1hrPriorMin' };
AUCTable(1,1:13) = ControlVarNames;
OPTROCPTTable = AUCTable;

[trainInd,valInd,testInd] = dividerand(length(GoodRowsindices),.7,0.3,0); % 0.7, 0.3
GoldTraining = NewDataVols(trainInd,2)+1;  % it needs to have positive categories: 0->1; 1->2;
GoldValid =  NewDataVols(valInd,2)+1;

[ ControlVars] = ConvertControlVarsToMat(AllRows,data);
ControlVarsGood = ControlVars(GoodRowsindices,2:end) ;% not dealing with ICUstay now
ControlVarsGood = [ControlVarsGood NaN*ones(length(GoodRowsindices), 3)];

for indTime = 1: length(TimeThresholdRange)% 1: length(TimeThresholdRange)
    lindx = sub2ind([length(AvgVolumeThresholdRange),length(TimeThresholdRange)],1,indTime); % row, column same order as size
    
    TimeDuration = TimeThresholdRange(indTime); AUCTable{indTime+1,1} = TimeDuration; OPTROCPTTable{indTime+1,1} = TimeDuration;
    
    % adding on to the control variables matrix the time dependent control
    % vars - this would change every round.
    
    ControlVarsGood(:,9) = totalFluidPriorMin(:,2+ lindx);
    ControlVarsGood(:,10) = TFVasopressorPriorBelowMin(:,2+ lindx);
    ControlVarsGood(:,11) = avgMAP1hrPriorBelowMin(:,2+ lindx);
    
    % plain model with no control variables
    
    [B] = mnrfit(NewDataVols(trainInd,2+ lindx),GoldTraining,'model','ordinal'); % B = model, coefficients including a constant
    [pihat] = mnrval(B,NewDataVols(valInd,2+ lindx),'model','ordinal');% predicted probabilities for each multinomial category
    
    [X,Y,T,AUC] = perfcurve(GoldValid,pihat(:,2),[2]); %% Probabilty or True, T = thresholds, X, Y : 1- spec , sens?
    
    [distJpoint,minInd] = nanmin(sqrt((X-0).^2  + (Y-1).^2));
   
    % Get CI for AUC and distance from J point
    
    [CIAUC_cont, CIdistJpoint_cont] = Compute95CI_AUC_distJPoint(GoldValid,pihat(:,2),200, 5); %%% change back to 200
    AUCTable{indTime+1,2}  = ['(',num2str(CIAUC_cont(1)),', ',num2str(AUC),', ',num2str(CIAUC_cont(2)),')'];
    OPTROCPTTable{indTime+1,2} = ['(',num2str(CIdistJpoint_cont(1)),', ',num2str(distJpoint),', ',num2str(CIdistJpoint_cont(2)),')'];
    
    figure();
    plot(X,Y, 'o');
    xlabel('False positive rate'); ylabel('True positive rate')
    title(['ROC for classification by logistic regression- No controls time duration ',num2str(TimeDuration)]);
    text([X(1) X(end)],[Y(1) Y(end)],{num2str(T(1)),num2str(T(end))}) % label the first last and best cutoff
    text(X(minInd),Y(minInd),num2str(T(minInd)));
    print(['ROC plots - No controls time duration ',num2str(TimeDuration)],'-dpng','-r300');
    XNoCtrl = X; YNoCtrl = Y;
    
    % making new folder
    newFolderName = ['dur',num2str(TimeDuration)];
    dirnew = mkdir(newFolderName);
   
    
    for controlInd = 1:11 
        
        [B] = mnrfit([NewDataVols(trainInd,2+ lindx) ControlVarsGood(trainInd,controlInd)],GoldTraining,'model','ordinal');
        [pihat] = mnrval(B,[NewDataVols(valInd,2+ lindx) ControlVarsGood(valInd,controlInd)],'model','ordinal');
        
        [X,Y,T,AUC] = perfcurve(GoldValid,pihat(:,2),[2]); %% Probabilty or True, T = thresholds, X, Y : 1- spec , sens?
       
        [distJpoint,minInd] = nanmin(sqrt((X-0).^2  + (Y-1).^2));
        % Get CI for AUC and distance from J point
        [CIAUC, CIdistJpoint] = Compute95CI_AUC_distJPoint(GoldValid,pihat(:,2),200, 5);
        
        % See if there is any overlap with 'no control'
        
        overlapAUCTF = ((CIAUC_cont(1) > CIAUC(2)) |(CIAUC_cont(2) < CIAUC(1)));
        overdistJpointTF = ((CIdistJpoint_cont(1) > CIdistJpoint(2)) |(CIdistJpoint_cont(2) < CIdistJpoint(1)));
  
        
        AUCTable{indTime+1,2 + controlInd }  = ['(',num2str(CIAUC(1)),', ',num2str(AUC),', ',num2str(CIAUC(2)),', ',num2str(overlapAUCTF),')'];
        OPTROCPTTable{indTime+1,2 + controlInd }= ['(',num2str(CIdistJpoint(1)),', ',num2str(distJpoint),', ',num2str(CIdistJpoint(2)),', ',num2str(overdistJpointTF),')'];
        
         cd (newFolderName) ;
        % Plotting
        figure();
        plot(X,Y, 'o'); hold on; plot(XNoCtrl,YNoCtrl,'r-'); legend({ControlVarNames{2+controlInd},'No ctrl'});
        legend('Location', 'SouthEast');
        xlabel('False positive rate'); ylabel('True positive rate')
        title(['ROC for classification by logistic regression- ', ControlVarNames{2+controlInd}, ' time duration ',num2str(TimeDuration)]);
        text([X(1) X(end)],[Y(1) Y(end)],{num2str(T(1)),num2str(T(end))}) % label the first last and best cutoff
        text(X(minInd),Y(minInd),num2str(T(minInd)));
        print(['ROC plots -  ', ControlVarNames{2+controlInd}, ' time duration ',num2str(TimeDuration)],'-dpng','-r300');
        cd ..
        
    end
    
    close all;
end
save('ROCValues_FixedAKIUO_NotTrunc.mat','AUCTable','OPTROCPTTable');
%save('ROCValues_Trunc.mat','AUCTable','OPTROCPTTable');

% Closing pools
parpool('close');


%% Histograms of the start times of the min of the urine output. It is in seconds. Not updated!
lindx = sub2ind([length(AvgVolumeThresholdRange),length(TimeThresholdRange)], 8,1);
%TimeThresholdRange = [2,4,6,8,10,12,14,16,18,20,22,24];% colSub
%AvgVolumeThresholdRange = [0: 0.1: 1];%rowSub

figure();
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz);
hist(NewDataStartTimeMins(:,2 + lindx)/3600,100);
hold on
plot([48,48],[0, nanmax(NewDataStartTimeMins(:,2 + lindx)/3600)],'r')
set(gca,'XTick',48,'XTickLabel','48')
xlabel('Start time first time below theshold point 7 and 2 hours')
ylabel('Count')
print('Start time first time below theshold point 7 and 2 hours','-dpng','-r300');



%% For each ICU calculate the confusion matrix for the original thresholds
% 1: CSRU
% 2: SICU
% 3: CCU
% 4: MICU + FICU
% ICUNumber
% 3: Gold , 4: Calculated.
%% All
%AvgVolumeThresholdRange = [0: 0.1: 1];%rowSub
%TimeThresholdRange = [2,4,6,8,10,12,14,16,18,20,22,24];% colSub

lindx = sub2ind([length(AvgVolumeThresholdRange),length(TimeThresholdRange)],6,3);
OrigThresholdComboIndex = 2 + lindx;
'All'

[CM,order] = confusionmat(NewData(:,OrigThresholdComboIndex),NewData(:,2),'ORDER',[1,0])
ConfusionMatrixInfo(CM);
% Venn diagram
A = sum(CM,2); A = A(1)
B = sum(CM,1);B = B(1)


venn([B, A], [CM(1,1)])
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz);

print('Venn diagram','-dpng','-r300');

%% 1: CSRU
'CSRU'

ICUColumn = ConvertToCell(data,14,AllRows);

TFMatchesICU = strcmp('CSRU',ICUColumn);
RowsMatchingICU = find(TFMatchesICU==1);
[CommonRows,ia,ib] = intersect(RowsMatchingICU,GoodRowsindices);

[CM,order] = confusionmat(NewData(ib,OrigThresholdComboIndex),NewData(ib,2),'ORDER',[1,0])
ConfusionMatrixInfo(CM);


%% 2: SICU
'SICU'

TFMatchesICU = strcmp('SICU',ICUColumn);
RowsMatchingICU = find(TFMatchesICU==1);
[CommonRows,ia,ib] = intersect(RowsMatchingICU,GoodRowsindices);

[CM,order] = confusionmat(NewData(ib,OrigThresholdComboIndex),NewData(ib,2),'ORDER',[1,0])

ConfusionMatrixInfo(CM);



%% 3: CCU
'CCU'

TFMatchesICU = strcmp('CCU',ICUColumn);
RowsMatchingICU = find(TFMatchesICU==1);
[CommonRows,ia,ib] = intersect(RowsMatchingICU,GoodRowsindices);

[CM,order] = confusionmat(NewData(ib,OrigThresholdComboIndex),NewData(ib,2),'ORDER',[1,0])

ConfusionMatrixInfo(CM);

%% 4: MICU + FICU

' MICU + FICU'

TFMatchesICU = (strcmp('MICU',ICUColumn) | strcmp('FICU',ICUColumn));
RowsMatchingICU = find(TFMatchesICU==1);

[CommonRows,ia,ib] = intersect(RowsMatchingICU,GoodRowsindices);

[CM,order] = confusionmat(NewData(ib,OrigThresholdComboIndex),NewData(ib,2),'ORDER',[1,0])
ConfusionMatrixInfo(CM);

%% Plot the min urine output (mL/kg/hr) of six hour intervals
close all

lindx = sub2ind([length(AvgVolumeThresholdRange),length(TimeThresholdRange)],6,3)
OrigThresholdComboIndex = 2 + lindx;

figure()
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
hist(NewDataVols(:,OrigThresholdComboIndex),10000);
xlim([0,10])
xlabel('Min urine output (mL per kg per hr) of six hour intervals');
ylabel('Count')
hold on
line([0.5 0.5],[0,70],'Color','red');
text(0.5,-30,'0.5')

title('Histogram of min urine output (mL per kg per hr) of six hour intervals');
print('Histogram of min urine output (mL per kg per hr) of six hour intervals','-dpng','-r300');








