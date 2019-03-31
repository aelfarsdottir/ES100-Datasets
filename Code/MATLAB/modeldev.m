%% Aldís Elfarsdóttir
% February 15-26

%% Consolidating model dev scripts including import function for training validation csv's

%% start
close all;
clc;

%% load structures of training and validation data
cd '/Users/Aldis/Documents/MATLAB/ES100/Dataset4/Preprocessed/021319'
load workspace022919_1030pm.mat % shared via GOOGLE DRIVE due to size
% https://drive.google.com/open?id=1WQeSVsohE-RZ6Zgh-qnToohKr_5el82F

%% initialize switch variables
structset = FULL;       % the structure containing training and validation data
inputvars = '14var';    % a tag for the number of input variables in this iteration of model development
structname = 'full';    % string name for the structure
makefigures = 'yes';    % suppress figure generation for linear simulations when sweeping over multiple sets etc.

%% make_iddata.m (system identification object)

% response variable for prediction
yvar = {'roomC'}; % may change in future

% input variables (to be partitioned into disturbance/control)
uvar = {'positR','positS','heatvalve','slabC','humid','co2','occ', ...
    'windir', 'rain','oatC','temp31','temp33','temp23','temp22'};

% to iterate over atrain, avalid, btrain, bvalid training sets in the struct stored & loaded in the workspace essentials.mat file
cellnames = fieldnames(structset);

%% MAKE IDDATA OBJECT FOR SYSTEM IDENTIFICATION USING N4SID

for cn = 1 : length(cellnames) % atrain, avalid, btrain, bvalid, before, after, fullrow
    
    % use these corenames to differentiate the iddata sets
    coreFileName = cellnames{cn}; % after, before, btrain etc.
    
    % iddataobj to function as training and validation set
    set = structset.(coreFileName);                     % extract specific training set
    uvars = zeros(length(set.time),length(uvar));       % compile vectors of input variables
    for i = 1:length(uvar)
        uvars(:,i) = set.(uvar{i});
    end
    iddataobj.(coreFileName) = iddata(set.(yvar),uvars);        % Input and output data for N4SID
    iddataobj.(coreFileName).Ts = 1; % 1 second interpolation   % Time step of interpolation / prediction
    iddataobj.(coreFileName).TimeUnit = 'minutes';              % Units
    iddataobj.(coreFileName).InputName = uvar; % multiple input % String of input names
    iddataobj.(coreFileName).OutputName = yvar; % single output % String of output names
    
    % confirmation
    fprintf(1, 'iddataobj_trainvalid: %s\n', coreFileName);     % print an update at the command line
    
end

%% RUN N4SID ON IDDATA OBJECT

order = [1,2,3]; % model order

for o = 1:length(order)
    
    ord = strcat('order',num2str(order(o))); % struct location
    init = strcat(ord,'x0'); % initial condition
    
    for cn = 1 : length(cellnames)-1
        
        coreFileName = cellnames{cn}
        
        % n4sid state space approx from training data (all corefiles)
        [MODEL.(coreFileName).(ord),...
            MODEL.(coreFileName).(init)] = ...
            n4sid(iddataobj.(coreFileName), order(o));
        
    end
    
end

%% COLLECT FIT TO VALIDATION DATA

tset = {'atrain','avalid','btrain','bvalid','fullrow'};
pem = {'n4sid'};
order = {'order1','order2','order3'};
makefigures = 'no';

for ts = 1 : length(tset)              % sweep training sets
    for pm = 1:length(pem)             % pm for parameter estimation method
        for od = 1:length(order)       % sweep orders
            
            % update model to sim with correct type and order
            modeltosim = MODEL.(tset{ts}).(order{od}); % for compare & lsim
            modelinit = MODEL.(tset{ts}).(strcat(order{od},'x0')); % for linear sim
            
            % COMPARISON & PREDICTION HORIZONS
            % for both lsim and compare, run through test data sets
            
            for vs = 1:length(tset)
                
                % validation data
                validata = IDDATAOBJ.(tset{vs}); 
                
                % comparison prediction horizon
                steps = [1,5,15,60]; % because these are by second...
                modelfit = zeros(1,length(steps));
                for st = 1:length(steps)
                    [~,fitti,~] = compare(validata, modeltosim, steps(st));
                    modelfit(st) = fitti; % for filing the report away
                end
                PREDHORZ.(tset{ts}).(order{od}).(tset{vs}) = modelfit
                
                % linear simulation
                % time and inputs from the validation set
                setv = FULL.(tset{vs});
                tvect = setv.time;
                uvect = zeros(length(modeltosim.InputName), length(tvect));
                for n=1:length(modeltosim.InputName) % number of inputs
                    uvect(n,:) = setv.(modeltosim.InputName{n})';
                end
                
                [ysim, tsim, xsim] = lsim(modeltosim,uvect',1:length(tvect),modelinit); % can it be allowed to initialize itself with its own initial state value... of course. it's the inputs that will be different!
                resid = setv.roomC - ysim;
                totdiff = setv.roomC - mean(setv.roomC);
                SSE = sum(resid.^2)
                SST = sum(totdiff.^2) % https://en.wikipedia.org/wiki/Total_sum_of_squares
                Rsq_tv = 1-(SSE/SST)
                Rsq.(tset{ts}).(order{od}).(tset{vs}) = Rsq_tv;

            end
        end
    end
end

%% Evaluate model fit to estimation data, and compile into a results table
% numcombo = length(pem)*length(tset)*length(order)*length(timehorz);
rowname = {};
fitpercent = [];
MSE = [];
AIC = [];
tset = {'atrain','avalid','btrain','bvalid','fullrow'};
order = {'order1','order2','order3'};
pem = {'n4sid'};
for pm = 1:length(pem)              % sweep parameter estimation method
    for ts = 1:length(tset)         % sweep training set
        for od = 1:length(order)
            
            % the Report in the MODEL struct gives fitestim values
            REFTABLE = MODEL.(tset{ts}).(order{od}).Report.Fit;

            % extract column values to add to MODELDEVresults table
            rowname = vertcat(rowname,strcat(tset{ts},'-',pem{pm},'-',order{od}));
            fitpercent = vertcat(fitpercent,REFTABLE.FitPercent);
            MSE = vertcat(MSE,REFTABLE.MSE);
            AIC = vertcat(AIC,REFTABLE.AIC);
            
            
        end
    end
end

MODELDEVresults.fitestim = table(rowname, fitpercent, MSE, AIC);

%% Save fit validation results into table
% numcombo = length(pem)*length(tset)*length(order)*length(timehorz);
rowname = {};
r2atrain = [];
r2avalid = [];
r2btrain = [];
r2bvalid = [];
r2fullrow = [];
r2before = [];
r2after = [];
tset = {'atrain','avalid','btrain','bvalid','fullrow'};
order = {'order1','order2','order3'};
pem = {'n4sid'};
for pm = 1:length(pem)              % sweep parameter estimation method
    for ts = 1:length(tset)         % sweep training set
        for od = 1:length(order)
            
            % the Report in the MODEL struct gives fitestim values
            REFTABLE = Rsq.(tset{ts}).(order{od});

            % extract column values to add to MODELDEVresults table
            rowname = vertcat(rowname,strcat(tset{ts},'-',pem{pm},'-',order{od}));
            r2atrain = vertcat(r2atrain,REFTABLE.atrain);
            r2avalid = vertcat(r2avalid,REFTABLE.avalid);
            r2btrain = vertcat(r2btrain,REFTABLE.btrain);
            r2bvalid = vertcat(r2bvalid,REFTABLE.bvalid);
            r2after = vertcat(r2after,REFTABLE.after);
            r2before = vertcat(r2before,REFTABLE.before);
            r2fullrow = vertcat(r2fullrow,REFTABLE.fullrow);
        end
    end
end

MODELDEVresults.fitvalid = table(rowname, r2atrain, r2avalid, r2btrain, r2bvalid, r2after, r2before, r2fullrow);

%% From given folder, import files and save as structure .mat 
% only need to run this script once, and it saves all variables into .mat
myFolder = '/Users/Aldis/Documents/MATLAB/ES100/Dataset4/Preprocessed/021319';
full = importdset_loadfolder(myFolder);

function FULL = importdset_loadfolder(myFolder)
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
save dataset4_021319_fullrow.mat full % for later imports (as at beginning of this script)
end
