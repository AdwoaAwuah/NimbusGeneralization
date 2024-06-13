%% Adding Norm to adaptData

% Load data and Find norms for the entire time courses
%This code will find Euclidean norm at each step for the entire time courses
%Created by DMM0 5/10/2022

% 1) load subjects
% 2) EMG normalization of baseline
% 3) Computing Stride by stride norm
% 4) Compute bias removed stride by stride norm
% 5) Saving params file 

%clear;clc; close all
%% 1: load and prep data
subID1= 'NTS_01';
load([subID1, 'params.mat'])

subID = {[subID1, 'params.mat']}; %Change SubID to a cell array to avoid brace indexing error when creating a group

%
%% 2:  EMG normalization of baseline
group=adaptationData.createGroupAdaptData(subID); % create group adaptation data from the subject.


muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'TFL', 'GLU', 'HIP'};
% muscleOrder={'TA'};
n_muscles = length(muscleOrder); %the number of muscles to plot

 ep=defineEpochs_regressionYA('nanmean'); % the epochs of interest
%  ep=defineRegressorsNimbusFeedback('nanmean');
 
refEp= defineReferenceEpoch('TM base',ep); %define reference or baseline for the rest of the epochs

epochOfInterest={'Adaptation','Adaptation_{early}','OG base'}; % This line chooses the epochs that want to get data from

newLabelPrefix = defineMuscleList(muscleOrder); 

%adaptData = adaptData.normalizeToBaselineEpoch(newLabelPrefix,refEp); %normalizing the data

group = group.normalizeToBaselineEpoch(newLabelPrefix,refEp);%commented above when changed to "group" 

%ll=adaptData.data.getLabelsThatMatch('^Norm');

ll=group.adaptData{1}.data.getLabelsThatMatch('^Norm');%changed to this

l2=regexprep(regexprep(ll,'^Norm',''),'_s','s');

%adaptData=adaptData.renameParams(ll,l2);
group=group.renameParams(ll,l2);

newLabelPrefix = regexprep(newLabelPrefix,'_s','s');

%
%

%Plotting Checkerboards

% fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
% ph=tight_subplot(1,length(epochOfInterest),[.03 .005],.04,.04);
% 
% flip=1; %1 for individual leg analysis and 2 for asymmetric (R-L)
% if flip==1
%     n=2;
%     method='IndvLegs';
% else
%     n=1;
%     method='Asym';
% end
% 
% removeBias=1; %Flag for bias removal This can be turned on to remove the bias before model fitting 
% 
% 
% for l=1:length(epochOfInterest)
% 
%     ep2=defineReferenceEpoch(epochOfInterest{l},ep);
% 
% 
%   if removeBias == 1
% 
% 
%     [~,~,~,Datacheck{l}]=group.plotCheckerboards(newLabelPrefix,ep2,fh,ph(1,l),refEp,flip);
%      
%   end
% 
% end

% Removing bad muscles
   group= RemovingBadMuscleToSubj(group); %Removing bad muscles.

%% 2. Norm Stride by Stride

% Defining needed variables
data=[];
temp=[];
aux1=[];

Subj = group.adaptData{1}; %Dummy variable
%Subj = adaptData

for i = 2:numel(newLabelPrefix) %loop on the all the muscles
    
    DataIdx=find(contains(Subj.data.labels, {[newLabelPrefix{i}, ' ']})); %Find data index (row where the muscles are)
    
    if length(DataIdx)<12 % In case the code does not grab all the muscles
        %(It should be 12 gaits phases of the gait cycle)
        DataIdxlast=DataIdx(end)+[1:3];
        DataIdx= [DataIdx; DataIdxlast'];
    end
    
    
    data=[data Subj.data.Data(:,DataIdx(3:6))]; %Concatenating all the muscle data
    data(isnan(data))=0; % if nan set to zero the norm function cant work with nan
    
end

data(isnan(data))=0;
dataAsym=data-fftshift(data,2); % For asymmetry measure sustract the second part of the matrix
dataAsym=dataAsym(:,1:size(dataAsym,2)/2,:); % Getting only the difference between legs
temp(:,1)=vecnorm(data'); % getting the norm
temp(:,2)=vecnorm(dataAsym'); % getting norm asymmetry value

%         aux1=find(temp(:,1)>50);
%         temp(aux1,:)=nan;

adaptData.data=adaptData.data.appendData(temp,{'NormEMG','NormEMGasym'},...
     {'Norm of all the muscles','Norm asym of all the muscles'}); % Adding parameter for to adaptData
% 
% group.adaptData{1}=group.adaptData{1}.appendData(temp,{'NormEMG','NormEMGasym'},...
%     {'Norm of all the muscles','Norm asym of all the muscles'}); % Adding parameter for to adaptData

%% 2. Norm Stride by Stride with baseline remove 

 ep=defineEpochs_regressionYA('nanmean');
% ep=defineRegressorsNimbusFeedback('nanmean'); %Define epochs of interest 
refEpTR= defineReferenceEpoch('TM base',ep); % defining  Treadmill baseline 
refEpOG= defineReferenceEpoch('OG base',ep);  % defining  overgrounds baseline 

padWithNaNFlag=true; % Fill with Nan in case that we dont have enought strides

%[OGref]=adaptData.getPrefixedEpochData(newLabelPrefix,refEpOG,padWithNaNFlag); % getting overgound baseline data
[OGref]=group.adaptData{1}.getPrefixedEpochData(newLabelPrefix,refEpOG,padWithNaNFlag);
OGref=squeeze(OGref); 
OGrefasym=OGref-fftshift(OGref,1); % getting OG base for the asymmetry parameter 
OGrefasym=OGref(1:size(OGref,1)/2,:,:);

%[TRref]=adaptData.getPrefixedEpochData(newLabelPrefix,refEpTR,padWithNaNFlag); % getting treadmill baseline data
[TRref]=group.adaptData{1}.getPrefixedEpochData(newLabelPrefix,refEpTR,padWithNaNFlag); 

TRref=squeeze(TRref);
TRrefasym=TRref-fftshift(TRref,1);  % getting TM base for the asymmetry parameter 
TRrefasym=TRref(1:size(TRref,1)/2,:,:);

%Defining needed dummy variables 
data=[];
temp=[];
data3=[];
data3asym=[];

Subj = group.adaptData{1};


for i = 1:numel(newLabelPrefix) %loop on the all the muscles
    
    DataIdx=find(contains(Subj.data.labels, {[newLabelPrefix{i}, ' ']}));
    if length(DataIdx)<12
        DataIdxlast=DataIdx(end)+[1:3];
        DataIdx= [DataIdx; DataIdxlast'];
    end
    
    
    data=[data Subj.data.Data(:,DataIdx)];
    data(isnan(data))=0;
    
end

trial=find(contains(Subj.data.labels, {'trial'}));
tt=unique(Subj.data.Data(:,trial));

for t=1:length(tt) % loop on all the trials 
    
    zz=tt(t);
    aux2=[];
    aux3=[];
    
    if find(contains(Subj.data.trialTypes(zz), {'OG'} )) %IF they are type OG remove OG baseline 
        
        Idx = find(Subj.data.Data(:,trial)==zz);
        aux2=data(Idx,:)';
        data2= aux2-OGref(:,1);
        
        aux3=aux2-fftshift(aux2,1); % For asymmetry measure sustract the second part of the matrix
        aux3=aux3(1:size(aux3,1)/2,:,:);
        
        data2asym=aux3-OGrefasym(:,1);

        
    else  %If they are type TM remove TM baseline 
        
        Idx = find(Subj.data.Data(:,trial)==zz);
        aux2=data(Idx,:)';
        data2= aux2-TRref(:,1);
        
        aux3=aux2-fftshift(aux2,1);% For asymmetry measure subtract the second part of the matrix
        aux3=aux3(1:size(aux3,1)/2,:,:);
        
        data2asym=aux3-TRrefasym(:,1);
        
        
    end
    

    data3=[data3 data2];
    data3asym=[data3asym data2asym];
    
end

data3(isnan(data3))=0;
data3asym(isnan(data3asym))=0; %converting nans to 0
temp(:,1)=vecnorm(data3);
temp(:,2)=vecnorm(data3asym);
%         aux1=find(temp(:,1)>50);
%         temp(aux1,:)=nan;
aux1=adaptData.data.Data;
adaptData.data=adaptData.data.appendData(temp,{'UnBiasNormEMG','UnBiasNormEMGasym'},...
    {'Context specific unbais Norm of all the muscles','Context specific unbias Norm asym of all the muscles'});  % Adding parameter for to adaptData

%% Plot some of the parameters 

 params= {'UnBiasNormEMG'};
%  params= {'NormEMG'};
% adaptData.plotAvgTimeCourse(adaptData,params)

conditions={'OG base','Adaptation','Post 1'};
%  conditions={'Post 1'};
trialMarkerFlag=0; %1 if you want to separete the time course by trial, 0 to separate by condition 
indivFlag=0; %0 we are plotting one subject
indivSubs=0; %0 we are plotting one subject
colorOrder=[]; %Let the function take care of this at least you wanted in a specific set of color then by my guess and add the list here
biofeedback= 0; % At least that you are providing with biofeedback to the subject
removeBiasFlag=1; %if you want to remove bias 

figure;
adaptData.plotAvgTimeCourse(adaptData,params,conditions,5,trialMarkerFlag,indivFlag,indivSubs,colorOrder,biofeedback,removeBiasFlag)


hold on;


yline(0,LineWidth=2)
%group.adaptData{1}.plotAvgTimeCourse(group.adaptData{1},params,adaptData.metaData.conditionName,5)
%% Save params file 

   save([subID1 'paramsEMGnorm.mat'],'group','adaptData')
