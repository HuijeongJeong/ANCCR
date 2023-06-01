% simulate teleport task used in Kim et al., 2020 with ANCCR

clearvars; close all; clc;
rng(2)

% task parameter
ntrial = 1000;
ncue = [9,1];
cuerewdelay = [1,9]; % delay from last cue onset to reward
cuecuedelay = 1; % delay b/w cs1 and cs2
consumdelay = 6;
meanITI = [2,60]; % short and long ITI

% model parameter
samplingperiod = 0.2;
w = 0.5;
k = 1;
Tratio = 2;
alpha = 0.02;
alpha_r = 0.2;
threshold = 0.4;
minimumrate = 10^(-3);
beta = [zeros(1,9),0.5]';
maximumjitter = 0.1;

nExp = size(meanITI*2,1); % number of experiments

%%
darsp = cell(4,1);
nIter = 3;
for iiter = 1:nIter
    iiter
    for iITI = 1:2
        %% generate eventlog
        
        for itone = 1:2 % 1:increasing, 2:steady
            if iiter==1
                darsp{iITI+(itone-1)*2} = nan(nIter,max(ncue)+1);
            end
            
            if itone==1
                eventlog = simulateEventChain(ntrial, ncue(itone), nan, meanITI(iITI), ...
                    meanITI(iITI)*3, cuecuedelay, cuerewdelay(itone), 1, consumdelay);
            else
                eventlog = simulateEvents(ntrial, 1, 2, 1, nan, meanITI(iITI),...
                    meanITI(iITI)*3, cuerewdelay(itone), 1, consumdelay);
            end
            
            %% simulate
            % assume change in T
            IRI = cuecuedelay*(ncue(itone)-1)+cuerewdelay(itone)+meanITI(iITI)+consumdelay;
            [DA,ANCCR,PRC,SRC,NC,~,~,Mij,Mi] = calculateANCCR_v2(eventlog,IRI*Tratio,...
                alpha, k, samplingperiod,w,threshold,minimumrate,beta([1:ncue(itone),ncue(1)+1]),alpha_r,maximumjitter);
            incue = find(eventlog(:,1)==1);
            
            if itone==1
                darsp{iITI+(itone-1)*2}(iiter,:) = mean(DA(incue+[0:ncue(itone)]),1);
            else
                darsp{iITI+(itone-1)*2}(iiter,[1,end]) = mean(DA(incue+[0:ncue(itone)]),1);
            end
        end
    end
end

%%
ITItype = {'Short ITI', 'Long ITI'};
clr_light = {[0.6 0.6 0.6],[1 0.6 0.6]};
clr = {'k','r'};
for iITI = 1:2
    axes('Position',axpt(2,1,iITI,1));
    hold on;
for itone = 1:2
    plot(0:9,darsp{iITI+(itone-1)*2},'Color',clr_light{itone});
    errorbar(0:9,mean(darsp{iITI+(itone-1)*2},1),std(darsp{iITI+(itone-1)*2},[],1)/sqrt(nIter),clr{itone});
end
ylim([0 1.2])
xlim([-1 10])
title(ITItype{iITI})
plot([0 0 nan 9 9],[0 1.2 nan 0 1.2],'k:')
xlabel('Time from cue (s)')
if iITI==1
ylabel('simulated DA response')
end
end
%%

