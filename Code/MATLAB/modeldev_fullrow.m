%% Consolidating model dev and prelim control scripts
% Feb 15,2019 4:33PM-5:44PM
% goals are to:
    % 1) Consolidate script so functions are in one document -- DONE
    % 2) Model comparison plots
        % overlay models and calculate MSE on randomly generated roomC
        % validation
    % 3) Sensitivity analysis -- to do 2/24-2/25
        % add errors to input data and see if model can maintain MSE
    % 4) Comprehensive branch analysis -- IN PROGRESS 2/19
        % test how well shorter training models predict on validation- 
        % timespans (from fullrow) vs. longer training set models
            % we have avalid, atrain, bvalid, btrain (1/4 size)
            % before, after (1/2 size) 
            % fullrow (1 size), probably generate random verification set
        % test higher order models [C1 C2 C3 ...]*[x1 x2 x3 ...]' all added
        % up will give the y for each time step -- DONE 2/19
            % expand table of results to accomodate higher orders
        % test n4sid, ssest, cvx model dev... 
        % do variable selection tests on (each of 6) model and tabulate R^2
            % if R^2 goes down, don't remove the variable
            % if R^2 goes up, remove the variable
            % also try fitlm tests
        % generate random verification sets and track MSE in a table/plot

% Feb 24,2019 1:11PM
% with dataset5, goals are to:
    % 1) load the data
    % 2) plot results from experiments (break down by time interval and
    %    control method settings
    % 3) analyze results, note areas of improvement for control script and
    %    actual implementation
    % 4) add to written report, leaving room in tables for next round of
    %    results
    % 5) it's OK if it didn't work as expected
    
% Feb 26, 2019
% see if Lina can help fix the code for trial this Friday or Friday 15th
% lay out analysis for March 8 report

%% start
close all;
clc;
cd '/Users/Aldis/Documents/MATLAB/ES100/Dataset4/Preprocessed/021319'
load dataset4_021319_fullrow.mat
load workspace022919_1030pm.mat
% load workspace022010_1030pm.mat
%%
% initialize switch variables
structset = full;
inputvars = '14var';
structname = 'full';
makefigures = 'yes'; % suppress if not

%% make_iddata.m

% predicted variable
yvar = {'roomC'}; % may change in future

% input variables (to be partitioned into disturbance/control)
uvar = {'positR','positS','heatvalve','slabC','humid','co2','occ', ...
    'windir', 'rain','oatC','temp31','temp33','temp23','temp22'};

% what if we take out temp31?

cellnames = fieldnames(structset);
%% training and validation iddata objects
for o = 1 : length(cellnames)
    
    % use these corenames to differentiate the iddata sets
    coreFileName = cellnames{o}; % after, before, btrain etc.
    
    % iddataobj to function as training and validation set
    set = structset.(coreFileName);
    uvars = zeros(length(set.Time),length(uvar));
    for i = 1:length(uvar)
        uvars(:,i) = set.(uvar{i});
    end
    iddataobj.(coreFileName) = iddata(set.(yvar{1}),uvars);
    iddataobj.(coreFileName).Ts = 1; % 1 minute interpolation
    iddataobj.(coreFileName).TimeUnit = 'minutes';
    iddataobj.(coreFileName).InputName = uvar; % multiple input
    iddataobj.(coreFileName).OutputName = yvar; % single output
    
    % confirmation
    fprintf(1, 'iddataobj_trainvalid: %s\n', coreFileName);
    
end
%% run_iddata
% Run IDDATA object (assuming all are big enough)
% just do n4sid and ssest same time
% switched "1" to "order"

order = 1; % model order
ord = strcat('order',num2str(order)); % struct location
init = strcat(ord,'x0'); % initial condition

for o = 1 : length(cellnames)

    coreFileName = cellnames{o}

    % n4sid state space approx from training data (all corefiles)
    [model.n4sid.(coreFileName).(ord),...
        model.n4sid.(coreFileName).(init)] = ...
        n4sid(iddataobj.(coreFileName), order);
%     [model.vn4sid.(coreFileName).(ord),...
%         model.vn4sid.(coreFileName).(init)] = ...
%         n4sid(iddataobj.v.(coreFileName), order); % redundant because I
%         made training and validation sets from each corefilename
    
    % ssest
    [model.ssest.(coreFileName).(ord),...
        model.ssest.(coreFileName).(init)] = ...
        ssest(iddataobj.(coreFileName), order);
%     [model.vssest.(coreFileName).(ord),...
%         model.vssest.(coreFileName).(init)] = ...
%         ssest(iddataobj.v.(coreFileName), order);
end
% change tn4sid to n4sid later
% was double counting

%% evaluate model fit to estimation data
% can simplify later with modeldev_031819 in 
% '/Users/Aldis/Documents/MATLAB/ES100/Dataset5/Preprocessed/031319'
% fit to estimation data (add it to the growing struct for iddata results)
for o = 1 : length(cellnames)

    coreFileName = cellnames{o}
    
    % temporary variables
    report.n4sid.('order1') = model.tn4sid.(coreFileName).('order1').Report.Fit;
    report.n4sid.('order2') = model.tn4sid.(coreFileName).('order2').Report.Fit;
    report.n4sid.('order3') = model.tn4sid.(coreFileName).('order3').Report.Fit;
    
    report.ssest.('order1') = model.tssest.(coreFileName).('order1').Report.Fit;
    report.ssest.('order2') = model.tssest.(coreFileName).('order2').Report.Fit;
    report.ssest.('order3') = model.tssest.(coreFileName).('order3').Report.Fit;
    
    iddataobj.report.(coreFileName).fitestim = {'Fitestim',...
        'n4sid-ord1','n4sid-ord2','n4sid-ord3',...
        'ssest-ord1', 'ssest-ord2','ssest-ord3';
        'FitPercent', report.n4sid.('order1').FitPercent,...
                        report.n4sid.('order2').FitPercent,...
                        report.n4sid.('order3').FitPercent,...
                        report.ssest.('order1').FitPercent,...
                        report.ssest.('order2').FitPercent,...
                        report.ssest.('order3').FitPercent;
        'MSE',report.n4sid.('order1').MSE,...
                        report.n4sid.('order2').MSE,...
                        report.n4sid.('order3').MSE,...
                        report.ssest.('order1').MSE,...
                        report.ssest.('order2').MSE,...
                        report.ssest.('order3').MSE;
        'AIC',report.n4sid.('order1').AIC,...
                        report.n4sid.('order2').AIC,...
                        report.n4sid.('order3').AIC,...
                        report.ssest.('order1').AIC,...
                        report.ssest.('order2').AIC,...
                        report.ssest.('order3').AIC;
        }  % metrics for fit to estimation data
end
%% fit to validation data
% this is where we do 
for c = 1 : length(cellnames)

    makefigures = 'no';
    coreFileName = cellnames{c}
    
    % validate each core file on everything that it wasn't built on
    % eliminates the need for case switching! phew
    % the order can become a case as well, in a loop
    types = {'tn4sid','tssest'}; % change to n4sid, ssest after re-running script with just training iddataobjs
    for m = 1:length(types) % m for parameter estimation method
        
        orders = {'order1','order2','order3'};
        for o = 1:length(orders)
            %% update model to sim with correct type and order
            modeltosim = model.(types{m}).(coreFileName).(orders{o}); % for compare & lsim
            modelinit = model.(types{m}).(coreFileName).(strcat(orders{o},'x0')); % for linear sim
            
            %% LINEAR SIMULATION
            % pass in model data
            A = modeltosim.A;
            B = modeltosim.B;
            C = modeltosim.C;
            D = modeltosim.D;
            
            %% COMPARISON & PREDICTION HORIZONS
            datasets = {'before','btrain','bvalid','after','atrain','avalid','fullrow'};
            
            %% for both lsim and compare, run through valid data sets
            validset = datasets;
            indices = find(~cellfun(@isempty,validset')); % thank you MATLAB central: https://www.mathworks.com/matlabcentral/answers/42283-index-non-empty-cells-in-cell-array
            for i = 1:length(indices)
                validata = iddataobj.(validset{indices(i)}); % validation data
                
                % comparison prediction horizon
                steps = [1,5,15];
                modelfit = zeros(1,length(steps));
                for j = 1:length(steps)
                    [~,fitti,~] = compare(validata, modeltosim, steps(j));
                    modelfit(j) = fitti; % for filing the report away
                end
                predhorz.(types{m}).(coreFileName).(orders{o}).(validset{indices(i)}) = modelfit;
                
                % linear simulation
                % time and inputs from the validation set
                setv = structset.(validset{indices(i)});
                tvect = setv.Time;
                uvect = zeros(length(modeltosim.InputName), length(tvect));
                for n=1:length(modeltosim.InputName) % number of inputs
                    uvect(n,:) = setv.(modeltosim.InputName{n})';
                end
                
                [ysim, tsim, xsim] = lsim(modeltosim,uvect',1:1:length(tvect),modelinit); % can it be allowed to initialize itself with its own initial state value... of course. it's the inputs that will be different!
                resid = setv.roomC - ysim;
                totdiff = setv.roomC - mean(setv.roomC);
                SSE = sum(resid.^2);
                SST = sum(totdiff.^2); % https://en.wikipedia.org/wiki/Total_sum_of_squares
                Rsq_tv = 1-(SSE/SST)
                Rsq.(types{m}).(coreFileName).(orders{o}).(validset{indices(i)}) = Rsq_tv;

                if strcmp(makefigures,'yes')

                    cd '/Users/Aldis/Documents/MATLAB/ES100/Dataset4/Preprocessed/021319/Figures_021719'
                    
                    figure; % lsim
                    % lsim(modeltosim,uvect',1:1:length(tvect),modelinit)
                    % % gives too many plots
                    plot(tsim, ysim); hold on; % simulated data
                    plot(tsim, setv.roomC); % measured data
                    annotation('textbox',[.3 .3 .4 .5],'String',['Rsq = ', num2str(round(Rsq_tv,4))], 'FontSize', 16, 'FitBoxToText','on')
                    legend('simulated','measured')
                    ylabel('Room Temperature (ºC)')
                    title(['Linear Simulation: ',coreFileName,' vs. ', validset{indices(i)}])
                    filesave = strcat(structname, '_', inputvars, '_', types{m}, '_', orders{o}, '_', coreFileName,'_vs_',validset{indices(i)},'_lsim.png')
                    saveas(gcf,filesave)
                    
                    % resid
                    figure;
                    plot(tsim,resid)
                    legend('residuals')
                    title(['Residuals: ',coreFileName,' vs. ', validset{indices(i)}])
                    filesave = strcat(structname, '_', inputvars, '_', types{m}, '_', orders{o}, '_', coreFileName,'_vs_',validset{indices(i)},'_resid.png')
                    saveas(gcf,filesave)
                    
                    cd ..

                end
                % something very strange about plotting > 500 figures
                % 2 types, 3 orders, 7*7 setv pairs, 2 images --> 588 plots
                % stopped it at 538 or so because I thought it'd go to 700
            end
        end
    end  
end

%% into table
% each case has a different table!
for c = 1 : length(cellnames)

    coreFileName = cellnames{c};
    
    % IT'S THE BEGINNING BRACKET {} --> []
    iddataobj.report.(coreFileName).fitvalid = [...
        {'Train set', 'after',' ',' ',' ','atrain',' ',' ',' ',...
            'avalid',' ',' ',' ','before',' ',' ',' ','btrain',' ',...
            ' ',' ','bvalid',' ',' ',' ','fullrow',' ',' ',' '}; % had to convert to string
        'Metric:', repmat({'r2','1step','5step','15step'},1,7);
        ['n4sid-ord1',Rsq.tn4sid.(coreFileName).('order1').('after'),...
                        num2cell(predhorz.tn4sid.(coreFileName).('order1').('after')),...
                      Rsq.tn4sid.(coreFileName).('order1').('atrain'), ...
                        num2cell(predhorz.tn4sid.(coreFileName).('order1').('atrain')), ...
                      Rsq.tn4sid.(coreFileName).('order1').('avalid'),...
                        num2cell(predhorz.tn4sid.(coreFileName).('order1').('avalid')),...
                      Rsq.tn4sid.(coreFileName).('order1').('before'),...
                        num2cell(predhorz.tn4sid.(coreFileName).('order1').('before')),...
                      Rsq.tn4sid.(coreFileName).('order1').('btrain'), ...
                        num2cell(predhorz.tn4sid.(coreFileName).('order1').('btrain')),...
                      Rsq.tn4sid.(coreFileName).('order1').('bvalid'),...
                        num2cell(predhorz.tn4sid.(coreFileName).('order1').('bvalid')),...
                      Rsq.tn4sid.(coreFileName).('order1').('fullrow'),...
                        num2cell(predhorz.tn4sid.(coreFileName).('order1').('fullrow'))];
                        
        ['n4sid-ord2',Rsq.tn4sid.(coreFileName).('order2').('after'),...
                        num2cell(predhorz.tn4sid.(coreFileName).('order2').('after')),...
                      Rsq.tn4sid.(coreFileName).('order2').('atrain'), ...
                        num2cell(predhorz.tn4sid.(coreFileName).('order2').('atrain')), ...
                      Rsq.tn4sid.(coreFileName).('order2').('avalid'),...
                        num2cell(predhorz.tn4sid.(coreFileName).('order2').('avalid')),...
                      Rsq.tn4sid.(coreFileName).('order2').('before'),...
                        num2cell(predhorz.tn4sid.(coreFileName).('order2').('before')),...
                      Rsq.tn4sid.(coreFileName).('order2').('btrain'), ...
                        num2cell(predhorz.tn4sid.(coreFileName).('order2').('btrain')),...
                      Rsq.tn4sid.(coreFileName).('order2').('bvalid'),...
                        num2cell(predhorz.tn4sid.(coreFileName).('order2').('bvalid')),...
                      Rsq.tn4sid.(coreFileName).('order2').('fullrow'),...
                        num2cell(predhorz.tn4sid.(coreFileName).('order2').('fullrow'))];
                        
        ['n4sid-ord3',Rsq.tn4sid.(coreFileName).('order3').('after'),...
                        num2cell(predhorz.tn4sid.(coreFileName).('order3').('after')),...
                      Rsq.tn4sid.(coreFileName).('order3').('atrain'), ...
                        num2cell(predhorz.tn4sid.(coreFileName).('order3').('atrain')), ...
                      Rsq.tn4sid.(coreFileName).('order3').('avalid'),...
                        num2cell(predhorz.tn4sid.(coreFileName).('order3').('avalid')),...
                      Rsq.tn4sid.(coreFileName).('order3').('before'),...
                        num2cell(predhorz.tn4sid.(coreFileName).('order3').('before')),...
                      Rsq.tn4sid.(coreFileName).('order3').('btrain'), ...
                        num2cell(predhorz.tn4sid.(coreFileName).('order3').('btrain')),...
                      Rsq.tn4sid.(coreFileName).('order3').('bvalid'),...
                        num2cell(predhorz.tn4sid.(coreFileName).('order3').('bvalid')),...
                      Rsq.tn4sid.(coreFileName).('order3').('fullrow'),...
                        num2cell(predhorz.tn4sid.(coreFileName).('order3').('fullrow'))];
                        
        ['ssest-ord1',Rsq.tssest.(coreFileName).('order1').('after'),...
                        num2cell(predhorz.tssest.(coreFileName).('order1').('after')),...
                      Rsq.tssest.(coreFileName).('order1').('atrain'), ...
                        num2cell(predhorz.tssest.(coreFileName).('order1').('atrain')), ...
                      Rsq.tssest.(coreFileName).('order1').('avalid'),...
                        num2cell(predhorz.tssest.(coreFileName).('order1').('avalid')),...
                      Rsq.tssest.(coreFileName).('order1').('before'),...
                        num2cell(predhorz.tssest.(coreFileName).('order1').('before')),...
                      Rsq.tssest.(coreFileName).('order1').('btrain'), ...
                        num2cell(predhorz.tssest.(coreFileName).('order1').('btrain')),...
                      Rsq.tssest.(coreFileName).('order1').('bvalid'),...
                        num2cell(predhorz.tssest.(coreFileName).('order1').('bvalid')),...
                      Rsq.tssest.(coreFileName).('order1').('fullrow'),...
                        num2cell(predhorz.tssest.(coreFileName).('order1').('fullrow'))];
                        
        ['ssest-ord2',Rsq.tssest.(coreFileName).('order2').('after'),...
                        num2cell(predhorz.tssest.(coreFileName).('order2').('after')),...
                      Rsq.tssest.(coreFileName).('order2').('atrain'), ...
                        num2cell(predhorz.tssest.(coreFileName).('order2').('atrain')), ...
                      Rsq.tssest.(coreFileName).('order2').('avalid'),...
                        num2cell(predhorz.tssest.(coreFileName).('order2').('avalid')),...
                      Rsq.tssest.(coreFileName).('order2').('before'),...
                        num2cell(predhorz.tssest.(coreFileName).('order2').('before')),...
                      Rsq.tssest.(coreFileName).('order2').('btrain'), ...
                        num2cell(predhorz.tssest.(coreFileName).('order2').('btrain')),...
                      Rsq.tssest.(coreFileName).('order2').('bvalid'),...
                        num2cell(predhorz.tssest.(coreFileName).('order2').('bvalid')),...
                      Rsq.tssest.(coreFileName).('order2').('fullrow'),...
                        num2cell(predhorz.tssest.(coreFileName).('order2').('fullrow'))];
                        
        ['ssest-ord3',Rsq.tssest.(coreFileName).('order3').('after'),...
                        num2cell(predhorz.tssest.(coreFileName).('order3').('after')),...
                      Rsq.tssest.(coreFileName).('order3').('atrain'), ...
                        num2cell(predhorz.tssest.(coreFileName).('order3').('atrain')), ...
                      Rsq.tssest.(coreFileName).('order3').('avalid'),...
                        num2cell(predhorz.tssest.(coreFileName).('order3').('avalid')),...
                      Rsq.tssest.(coreFileName).('order3').('before'),...
                        num2cell(predhorz.tssest.(coreFileName).('order3').('before')),...
                      Rsq.tssest.(coreFileName).('order3').('btrain'), ...
                        num2cell(predhorz.tssest.(coreFileName).('order3').('btrain')),...
                      Rsq.tssest.(coreFileName).('order3').('bvalid'),...
                        num2cell(predhorz.tssest.(coreFileName).('order3').('bvalid')),...
                      Rsq.tssest.(coreFileName).('order3').('fullrow'),...
                        num2cell(predhorz.tssest.(coreFileName).('order3').('fullrow'))];
       ];
end

%% MODEL COEFFICIENTS

%% RESULTS - combined analytics table 3/3
% for results see iddataobj.report.(coreFileName) extract fitestim (transpose) 
% and fitvalid

bfite = iddataobj.report.before.fitestim';
btfite = iddataobj.report.btrain.fitestim';
bvfite = iddataobj.report.bvalid.fitestim';
afite = iddataobj.report.after.fitestim';
atfite = iddataobj.report.atrain.fitestim';
avfite = iddataobj.report.avalid.fitestim';

RESULTS.fitestim = [...
    [['Model'; repmat({'fullrow'},6,1)] iddataobj.report.fullrow.fitestim']; ...
    [repmat({'before'},6,1) bfite(2:end,:)]; ...
    [repmat({'btrain'},6,1) btfite(2:end,:)]; ...
    [repmat({'bvalid'},6,1) bvfite(2:end,:)]; ...
    [repmat({'after'},6,1) afite(2:end,:)]; ... 
    [repmat({'atrain'},6,1) atfite(2:end,:)]; ...
    [repmat({'avalid'},6,1) avfite(2:end,:)]];

writetable(cell2table(RESULTS.fitestim),'14var_results_fitestim_030319.csv')
% validation fit table

bfitv = iddataobj.report.before.fitvalid;
btfitv = iddataobj.report.btrain.fitvalid;
bvfitv = iddataobj.report.bvalid.fitvalid;
afitv = iddataobj.report.after.fitvalid;
atfitv = iddataobj.report.atrain.fitvalid;
avfitv = iddataobj.report.avalid.fitvalid;

RESULTS.fitvalid = [...
    [['Model'; 'Model'; repmat({'fullrow'},6,1)] iddataobj.report.fullrow.fitvalid]; ...
    [repmat({'before'},6,1) bfitv(3:end,:)]; ...
    [repmat({'btrain'},6,1) btfitv(3:end,:)]; ...
    [repmat({'bvalid'},6,1) bvfitv(3:end,:)]; ...
    [repmat({'after'},6,1) afitv(3:end,:)]; ... 
    [repmat({'atrain'},6,1) atfitv(3:end,:)]; ...
    [repmat({'avalid'},6,1) avfitv(3:end,:)]];

writetable(cell2table(RESULTS.fitvalid),'14var_results_fitvalid_030319.csv')

%% 2/20/19 reviewing ES155 A-BK controller
% 1st place for validation fits is fullrow ssest-ord2 and n4sid-ord3
% 2nd place is atrain with ssest-ord2 (91% fullrow, 80% before)

model1 = model.tssest.fullrow.order2
model2 = model.tssest.atrain.order2

A = model2.A;
B = model2.B;
C = model2.C;
D = model2.D;

eig = [-1+i, -1-i];
K = place(A,B,eig)
kr = ones(1,14)/(-C*inv(A-B*K)*B)

sys_ctrl = ss((A-B*K),B*kr,C,D);
sys = ss(A,B,C,D);
figure;
step(sys_ctrl)
legend('with control')

figure;
step(model2)
legend('without')
% occupancy is 16ºC jump with a step in occupancy... by 2000 minutes
% seems overall unrealistic
% what to do with this K?

%% control trial allow inputs

%% tracking error over all sensor time is probably too low for typical day
otrek_before = full.before.('roomC') - full.before.('setC');
otrek_btrain = full.btrain.('roomC') - full.btrain.('setC');
otrek_bvalid = full.bvalid.('roomC') - full.bvalid.('setC');
otrek_after = full.after.('roomC') - full.after.('setC');
otrek_atrain = full.atrain.('roomC') - full.atrain.('setC');
otrek_avalid = full.avalid.('roomC') - full.avalid.('setC');
otrek_fullrow = full.fullrow.('roomC') - full.fullrow.('setC');
avg_otrek = {'before',mean(otrek_before); 'btrain', mean(otrek_btrain); 'bvalid',mean(otrek_bvalid); ...
    'after', mean(otrek_after); 'atrain', mean(otrek_atrain); 'avalid', mean(otrek_avalid); ...
    'fullrow', mean(otrek_fullrow)} 
mean([avg_otrek{:,2}]) % about 0.58ºC off full data (omitting Oct2-26)

%% tracking error on a given 2 hour period (at beginning of dset)

trackingerrs = zeros(1,length(cellnames));

for c = 1 : length(cellnames)

    coreFileName = cellnames{c};
    
    segment = 1:120:size(full.(coreFileName).Time);

    for i = 1:length(segment)-1
        timeint = segment(i):segment(i+1);
        otrek.(coreFileName)(i) = mean(full.(coreFileName).('roomC')(timeint) - full.(coreFileName).('setC')(timeint));
    end
    
    coreFileName
    trackingerrs(c) = mean(otrek.(coreFileName))
    
end
% this only increases the mean mean tracking error to 0.5871



%% use switch function for input at command line
% model1 = model.tssest.fullrow.order2;
% model2 = model.tssest.atrain.order2;
% 
% ctrlmodel = model2;
% 
% A = ctrlmodel.A;
% B = ctrlmodel.B;
% Bctrl = B(:,1); % positR, positS
% F = B(:,2:end); % heatvalve, slabC, humid, co2, occ, windir, rain, oatC, temp31, temp33, temp23, temp22
% C = ctrlmodel.C;
% D = ctrlmodel.D;
% 
% ydes = 17.5; % 16ºC
% ctrlnames = ctrlmodel.InputName(1); % interested in controlling heatvalve and positR (didn't include positS yet)
% 
% yobs0 = 18.84; % FORGOT TO CHANGE THIS UNTIL 12:18
% xobs0 = C\yobs0; % essentially x(0) = y(0)/C
% % uobs0 = [2;2]; % cvx won't accept initial value for u
% dobs0 = [0	0	21.1	20.48	566.72	0	66	0	9.63	18.84	19.38	21.54	21.54]

% 1:40
% 0	0	21.1	20.48	566.72	0	66	0	9.63	18.84	19.38	21.54	21.54
%1:35
% 0	100	21.1	22.5	572.8	0	169	0	9.57	18.78	19.4	21.64	21.64
%1:30
% 0	100	21	22.5	574.72	0	304	0	9.55	18.88	19.36	21.6	21.6
%1:25
% 0	100	21	22.5	574.72	0	342	0	9.54	18.78	19.36	21.58	21.58
%1:20
% 0	100	21	22.42	584.96	0	333	0	9.51	18.98	1936	21.58	21.58
% 1:10
% 0	100	21	24.46	586.86	0	23	0	9.57	19.1	19.28	21.6	21.6
%1:00
% [0	100	20.9	22.7	580.8	0	323	0	9.69	19.08	19.34	21.38	21.38]
% 12:32
% [0	100	20.9	24.84	596.8	0	340	0	9.44	18.64	19.48	21.44	21.44]
% 12:25
% [0	100	20.9	24.8	608.96	0	328	0	9.29	18.74	19.5	21.56	21.56]
% 12:15
% 0	100	20.9	24.82	654.72	0	14	0	9.14	18.94	19.52	21.44	21.44
% 12:10
% [0	100	20.9	24.74	708.48	-1	16	0	8.8	19.04	19.66	21.46	21.46]
% 12:05
% [0	100	20.9	25.64	683	0	30	0	8.83	19.4	19.8	21.46	21.46]
% 12:00
% 11:47 [0 100	20.8	25.74	672.64	0	16	0	9.04	19.34	19.74	21.4	21.4]

% cpu trails: [100,23.5,64.04, 1.8586e3,-1, 359,0, 10.28, 20.84, 20.4622,19.6947,21.1021]; % update every 5 minutes
    % going through finding max(full.atrain.___)
    % heatvalve, slabC, humid, co2 (1.8586e3), occ, 
    % windir, rain, oatC, (mean(oatC) = 10.28, min = 0.83, max = 21.84)
    % temp31, temp33, temp23, temp22 (means)

    % zeros(1,12) % will transpose later - lowers tracking error likely bc
    % model is built on a lot of data that has a lot of zeros in
    % disturbances
 % why did I have a for loop here?
% uobs0 = [2;2]; % cvx won't accept initial value for u
% reformatted for screenshots for report
% model1 = model.tssest.fullrow.order2;

% MODEL USED ON FEB 22: model.tssest.atrain.order2
ctrlmodel = model.tssest.atrain.order2; % can choose different model later
% model.tn4sid.after.order2
% model.tn4sid.atrain.order2; % hey, looks better than ssest!
ord = 2

A = ctrlmodel.A;
B = ctrlmodel.B;
Bctrl = B(:,1); % positR
F = B(:,2:end); % heatvalve, positS, slabC, humid, co2, occ, windir, 
                % rain, oatC, temp31, temp33, temp23, temp22
C = ctrlmodel.C;
D = ctrlmodel.D;

ydes = 20; % was 20ºC
ctrlnames = ctrlmodel.InputName(1); % controlling positR

yobs0 = 18.84; % FORGOT TO CHANGE THIS UNTIL 12:18
xobs0 = C\yobs0; % essentially x(0) = y(0)/C
dobs0 = [0,0,21.1,20.48,566.72,0,66,0,9.63,18.84,19.38,21.54,21.54];

    timeahead = 5 % edit
    cvx_begin %quiet
    
    variable u(length(ctrlnames),timeahead) % 1x5 = [positR, positS]
    expression x(ord,timeahead) % 2 states, initializes a zero vector
    expression d(13,timeahead) % 14 vars - 2 control vars
    expression room(1,timeahead)
    x(:,1) = xobs0; % manual input for control trial
    d(:,1) = dobs0'; % assume same disturbance for 5 minutes
    
    for k = 1:timeahead-1 % iterate over time
        d(:,k+1) = d(:,k); % assume same disturbance over 5 minutes
        u(:,k) <= 100; % bounds on window opening extent
        u(:,k) >= 0; % lower bound
        
%         x(:,k) >= -20
%         x(:,k) <= 50 % if huge constraint, it's feasible of course!
        % how do you detect redundant constraints
        x(:,k+1) = A*x(:,k) + Bctrl*u(:,k) + F*d(:,k); % state space model
        
        room(k) = C*x(:,k+1)
    end
        room(end) = 10;


    % changing this metric of tracking error changes control output
%     y = C*x;
%     y1 = C*x(:,2); 
    
%     heatloss = f(u) % now you will have R term
    
    minimize norm(room - ydes) %+ heatloss % 11.48ºC is high
    cvx_end
    %%
    
% overall tracking error afterward
norm(y - ydes,2) % meaning ydes is higher than y1...
% it really shouldn't be undershooting the temperature if it understands
% how slowly the dynamics of the room change with time....
% round(mean(y - ones(1,length(timeahead))*ydes),2)
roofposit = round(u,1)
% mean(C*x - ydes) % 0.028ºC (but that's theoretical 90% improvement...)
% (C*x - ydes*ones(1,timeahead))

%% ALL DATA PLOT 3/7/19

% positR
% heatvalve, positS, slabC, humid, co2, occ, windir, 
% rain, oatC, temp31, temp33, temp23, temp22
dset = full.fullrow;
figure;
plot(dset.Time, dset.positR,'.', 'MarkerSize',12); hold on;
plot(dset.Time, dset.positS,'.', 'MarkerSize',12)
plot(dset.Time, dset.heatvalve,'.', 'MarkerSize',12)
plot(dset.Time, dset.slabC,'.', 'MarkerSize',12)
plot(dset.Time, dset.humid,'.', 'MarkerSize',12)
plot(dset.Time, dset.co2,'.', 'MarkerSize',12)
plot(dset.Time, dset.occ,'.', 'MarkerSize',12)
plot(dset.Time, dset.windir,'.', 'MarkerSize',12)
plot(dset.Time, dset.rain,'.', 'MarkerSize',12)
plot(dset.Time, dset.oatC,'.', 'MarkerSize',12)
plot(dset.Time, dset.temp31,'.', 'MarkerSize',12)
plot(dset.Time, dset.temp33,'.', 'MarkerSize',12)
plot(dset.Time, dset.temp23,'.', 'MarkerSize',12)
plot(dset.Time, dset.temp22,'.', 'MarkerSize',12)
ax = gca;
ax.FontSize = 15;
ylabel('Various units (see legend)')
legend({'Skylight (% open)','South (% open)','Heatvalve (% of time on)',...
    'Slab ºC','Humidity (%)','CO_2 (ppm)','Occupancy (-1 = occupied)',...
    'Wind direction (º)','Rain (-1 = rain)','Outdoor (ºC)','Lounge 31 ºC', ...
    'Office 33 (ºC)','Office 23 (ºC)','Office 22 (ºC)'},'location','eastoutside')
legend('boxoff')
ax.XTick = linspace(dset.Time(1), dset.Time(end), 4);
datetick('x','mm/dd', 'keepticks')
filesave = strcat('alldata_0828_1203.png')
saveas(gcf,filesave)

%% total combinations for a table of results (and plots to produce)
combos = zeros(1,30);
for i = 1:1:30
    combos(i) = nchoosek(30,i);
end
numvarcombos = sum(combos); % the possible variable combinations to select
    % how many rows of results if we have:
        % 6 models
        % 2 parameter estimation methods
        % up to 10 model orders
        % nchoosek 30 choose 30!... variables
        % 5-7 validation sets per test
        % 1-2 metrics to compare (R^2, MSE)
        % over 3-4 prediction horizon steps
numresultrows = 6*2*10*sum(combos)*7*2*4; % 7.22e12 rows of results

% number of plots to produce can start dividing
% if we do plots by model order, we can have 7.22e11 plots
% moral of story this is a lot of analysis... how long would it take??

%% From given folder, import files and save as mat 
% only need to run this script once, and it saves all variables into .mat
myFolder = '/Users/Aldis/Documents/MATLAB/ES100/Dataset4/Preprocessed/021319';
full = importdset_loadfolder(myFolder);

function full = importdset_loadfolder(myFolder)
% Check that folder exists.  
if ~isfolder(myFolder)
  errorMessage = sprintf('Error: Try no C:/. The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage)); % Warn if it doesn't.
  return;
end

% Import all csv files

filePattern = fullfile(myFolder, '*.csv'); % all csv files in folder
procfiles = dir(filePattern); % processed files

for p = 1 : length(procfiles) % before
    
    FileName = procfiles(p).name;
    coreFileName1 = strrep(string(FileName),'fullrow-',''); % extract the * from '*.csv' (me)
    coreFileName = strrep(string(coreFileName1),'.csv',''); % extract the * from '*.csv' (me)

    % read as table, save into struct, replace Var1 with 'Time'
    full.(coreFileName) = readtable(FileName);
    full.(coreFileName).Properties.VariableNames(1) = {'Time'};

    % Confirmation
    fprintf(1, 'table made for: %s\n', coreFileName);
end

% SAVE AS MATFILE
cd '/Users/Aldis/Documents/MATLAB/ES100/Dataset4/Preprocessed/021319'
save dataset4_021319_fullrow.mat full % for later imports
end


% taken out of predhorz comparisons:
% if we want to eliminate the ornes that we don't care about
% e.g. at least validset = erase(datasets,coreFileName);
% but good to compare to all datasets for completeness sake
%         switch coreFileName
%             case 'before'
%                 validset = erase(datasets,{'before','btrain','bvalid'});
%                 % regexp(datasets,'b*') gives indices of the b- ones
%             case 'btrain'
%                 validset = erase(datasets,{'before','btrain'});
%                 % because 'btrain' is a subset of 'before'
%             case 'bvalid'
%                 validset = erase(datasets,{'before','bvalid'});
%             case 'after'
%                 validset = erase(datasets,{'after','atrain','avalid'});
%             case 'atrain'
%                 validset = erase(datasets,{'after','atrain'});
%                 % because 'btrain' is a subset of 'before'
%             case 'avalid'
%                 validset = erase(datasets,{'after','avalid'});
%             case 'fullrow'
%                 validset = erase(datasets,coreFileName);
%         end