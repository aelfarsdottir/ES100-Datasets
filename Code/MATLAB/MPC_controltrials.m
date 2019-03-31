%% Aldís Elfarsdóttir
% March 13, 2019 

%% Script for model predictive control trials with fullrow-2nd order and atrain-2nd order
% Tested with Dr. Yang Zheng that YALMIP, CVX working

%% start
close all;
clc;

%% CONTROL MODEL
cd '/Users/Aldis/Documents/MATLAB/ES100/Dataset4/Preprocessed/021319'
load workspace_031319_essentials.mat % shared via DRIVE due to file size
% https://drive.google.com/open?id=1WQeSVsohE-RZ6Zgh-qnToohKr_5el82F

% switch depending on control trial 
pem = {'tn4sid'}; % t for training
tset = {'atrain'}; % can sweep over {'after','atrain','avalid','before','btrain','bvalid','fullrow'};
order = {'order2'}; % can sweep over {'order1','order2','order3'};
timehorz = 10; % can sweep over [1,2,5,10,15]


trialnum = 13;                          % UPDATE DURING CONTROL TRIAL
for th = 1:length(timehorz)             % can sweep time horizon
    for pm = 1:length(pem)              % can sweep parameter estimation method
        for ts = 1:length(tset)         % can sweep training set
            for od = 1:length(order)    % can sweep order
                
                % select model
                ctrlmodel = MODEL.(pem{pm}).(tset{ts}).(order{od}); % can choose different model later
                
                
                % coefficients
                A = ctrlmodel.A;
                B = ctrlmodel.B;    % partitioned as follows:
                ctrlnames = ctrlmodel.InputName(1); % controlling positR
                Bctrl = B(:,1);     % positR
                F = B(:,2:end);     % heatvalve, positS, slabC, humid, co2, occ, windir,
                                    % rain, oatC, temp31, temp33, temp23, temp22
                C = ctrlmodel.C;
                D = ctrlmodel.D;
                K = ctrlmodel.K;
                
                % time horizon
                timeahead = timehorz(th) % can sweep
                timehorizon = strcat('step',num2str(timehorz(th))); % indexing
                
                % desired temperature (setpoint)
                ydes = 20; % was 17.5ºC during control trial to satisfy oatC < setC < roomC
                
                % UPDATE DURING CONTROL TRIAL!
                yobs0 = 19.54:-.01:18.54; % make a guess for what the observed temperature will be. 
                % ^^ this will use actual sensor data when implemented in real-time interface!            
                yobs = yobs0(1:timeahead+1);
                dobs0 = [0	0	21.9	41.2	540.8	0	222	0	19.53	21.48	22.42	21.82	21.82]; # UPDATE EACH TIME
                
                % switch order and state initiation
                switch order{od}
                    case 'order1'
                        ord = 1;
                        xobs0 = yobs(1)./C(1); % for now
                    case 'order2'
                        ord = 2;
                        xobs0 = [yobs(1)./C(1); 0]; % essentially x(0) = y(0)/C
                    case 'order3'
                        ord = 3;
                        xobs0 = [yobs(1)./C(1);0;0]; % for now
                end
                
                
                % YALMIP syntax with Yang, CVX syntax + Ke(t) and state-initialization on own
                
                % variables
                u = sdpvar(1,timeahead);
                x = sdpvar(ord,timeahead+1);
                y = sdpvar(1,timeahead+1);
                eobs = sdpvar(1,timeahead+1);
                
                % constraints
                constr = []; % initialize constraint vector
                constr = [constr, x(:,1) == xobs0]; % + Bctrl*u(1)]; % + F*dobs0'];
                constr = [constr, y(1) == C*x(:,1)]; % not exactly equal to yobs(1)
                constr = [constr, eobs(1) == yobs(1) - y(1)]; % initiate pred-error term
                
                % step forward prediction using the state-space model
                for k = 2:timeahead+1
                    
                    % included K and error term (ypredicted - yobserved)
                    constr = [constr, x(:,k) == A*x(:,k-1) + Bctrl*u(k-1) + F*dobs0' + K*eobs(k-1)];
                    
                    % response variable
                    constr = [constr, y(k) == C*x(:,k) + eobs(k)];
                    
                    % definition of prediction error
                    constr = [constr, eobs(k) == yobs(k) - y(k)];
                    
                end
                
                constr = [constr, u <= 100*ones(1,timeahead), u >= zeros(1,timeahead)]; % realistic window opening
                constr = [constr, y(end) == ydes]; % end constraint
                
                % cost function
                cost = norm(y - ydes*ones(1,timeahead+1)); % 2-norm of tracking error (roomC - setC)
                
                % call CVX solvers
                optimize(constr,cost);
                
                % report results
                solu = value(u)
                soly = value(y);
                solx = value(x);
                sole = value(eobs);
                solc = value(cost)
                
%                 % when sweeping over parameter estimation method, training set, time horizon, etc....
%                 OPT.(pem{pm}).(tset{ts}).(order{od}).(timehorizon).solu = value(u);
%                 OPT.(pem{pm}).(tset{ts}).(order{od}).(timehorizon).soly = value(y);
%                 OPT.(pem{pm}).(tset{ts}).(order{od}).(timehorizon).solx = value(x);
%                 OPT.(pem{pm}).(tset{ts}).(order{od}).(timehorizon).sole = value(eobs);
%                 OPT.(pem{pm}).(tset{ts}).(order{od}).(timehorizon).solc = value(cost);
                
                cd '/Users/Aldis/Documents/MATLAB/ES100/Dataset4/Preprocessed/021319/Controltrial_031519_atrain'
                figure;
                % build plotting vectors
                plotx = zeros(1,timeahead*2 + 2); % +2 for starting and end points to connect line segments of window position
                ploty = zeros(1,timeahead*2 + 2);
                addedx = zeros(1);
                addedy = zeros(1,2);
                for r = 1:timeahead
                    addx = repmat(r,1,2);
                    addy = repmat(solu(r),1,2);
                    addedx = horzcat(addedx, addx);
                    addedy = horzcat(addedy, addy);
                end
                addx = timeahead+1;
                plotx = horzcat(addedx,addx);
                ploty = addedy;
                plot(plotx', ploty'); hold on; %... [0,0, solu(1),solu(1),solu(2),solu(2),solu(3),solu(3),... solu(4),solu(4),solu(5),solu(5)]
                plot(0:timeahead+1, [yobs(1), soly])
                xlabel('Minutes')
                ylabel('Temp ºC and % Open')
                legend('PositR (% open)','MPC-predicted roomC','Location','Northwest')
                ax = gca;
                ax.FontSize = 12;
                
                % save new file for each segment of the control trial for later reference
                filesave = strcat(['controlsim_031519_',pem{pm},'_',tset{ts},'_',order{od},'_',num2str(timeahead),'min_trial',num2str(trialnum),'.png'])
                saveas(gcf,filesave)
                
            end
        end
    end
end

% print out solution (% opening of roof skylight window) for implementation by minute
solu
