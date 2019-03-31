%% Consolidating model dev and prelim control scripts
% March 13, 2019 testing with Yang Zheng that YALMIP, CVX working
% Found that model needs to be redeveloped for more significant parameters
% can I just test to see if more drastic A's and B's work?
%% start
close all;
clc;
cd '/Users/Aldis/Documents/MATLAB/ES100/Dataset4/Preprocessed/021319'
load workspace_031319_essentials.mat

%% TO DO

% COEFFICIENTS = ??
% SEE WHICH MODEL IS BEST AT CONTROLLING! (IN SIM AND IN LAB)
% MODIFY SCRIPT TO TAKE USER INPUTS FOR OBSERVED DISTURBANCES AND
% TEMPERATURES EACH MINUTE

%% CONTROL MODEL
% MODEL USED ON FEB 22: model.tssest.atrain.order2
cd '/Users/Aldis/Documents/MATLAB/ES100/Dataset4/Preprocessed/021319'
load workspace_031319_essentials.mat

pem = {'tn4sid'};

% switch depending on avalid (1st trial) atrain (2nd trial)
tset = {'atrain'}; %{'after','atrain','avalid','before','btrain','bvalid','fullrow'};
order = {'order2'}; %{'order1','order2','order3'};
timehorz = 10; % [1,2,5,10,15]; % sweep % done with 1,2,5,


trialnum = 13;
for th = 1:length(timehorz)         % sweep time horizon
    for pm = 1:length(pem)              % sweep parameter estimation method
        for ts = 1:length(tset)         % sweep training set
            for od = 1:length(order)    % sweep order
                
                % select model
                ctrlmodel = MODEL.(pem{pm}).(tset{ts}).(order{od}); % can choose different model later
                
                
                % coefficients
                A = ctrlmodel.A;
                B = ctrlmodel.B; % partitioned below
                ctrlnames = ctrlmodel.InputName(1); % controlling positR
                Bctrl = B(:,1);  % positR
                F = B(:,2:end);  % heatvalve, positS, slabC, humid, co2, occ, windir,
                % rain, oatC, temp31, temp33, temp23, temp22
                C = ctrlmodel.C;
                D = ctrlmodel.D;
                K = ctrlmodel.K;
                
                % time horizon
                timeahead = timehorz(th) % sweep
                timehorizon = strcat('step',num2str(timehorz(th))); % indexing
                
                % desired temperature
                ydes = 20; % was 20ºC
                
                % CHANGE TO INPUTS!
                yobs0 = 21.54:-.01:20.54; % slowly observing a temperature increase
                yobs = yobs0(1:timeahead+1);
                dobs0 = [0	0	21.9	41.2	540.8	0	222	0	19.53	21.48	22.42	21.82	21.82];
                
                % switch order and state initiation
                switch order{od}
                    case 'order1'
                        ord = 1;
                        xobs0 = yobs(1)./C(1); % for now
                    case 'order2'
                        ord = 2;
                        xobs0 = [yobs(1)./C(1); 0]; % yobs(1)./C(1) SWITCHING THIS REDUCED COST TO 7!! % essentially x(0) = y(0)/C
                    case 'order3'
                        ord = 3;
                        xobs0 = [yobs(1)./C(1);0;0]; % guessing for now
                end
                
                
                % Yalmip with Yang + Ke(t) and xinit on own
                
                % variables
                u = sdpvar(1,timeahead);
                x = sdpvar(ord,timeahead+1);
                y = sdpvar(1,timeahead+1);
                eobs = sdpvar(1,timeahead+1);
                
                % constraints
                constr = [];
                constr = [constr, x(:,1) == xobs0];% + Bctrl*u(1)]; % + F*dobs0'];
                constr = [constr, y(1) == C*x(:,1)]; % not exactly equal to yobs(1)
                
                % y-yobs --> 41 cost, yobs-y --> 271 cost! switching xobs0 "arbitrary 0
                % assignment" reduced cost to ~7. yobs-y now gives slightly lower cost
                constr = [constr, eobs(1) == yobs(1) - y(1)]; % initiate pred-error term
                
                for k = 2:timeahead+1
                    
                    % included K and error term (ypredicted - yobserved)
                    constr = [constr, x(:,k) == A*x(:,k-1) + Bctrl*u(k-1) + F*dobs0' + K*eobs(k-1)];
                    
                    constr = [constr, y(k) == C*x(:,k) + eobs(k)];
                    
                    % new: add error term predicted minus measured (let measured be yobs(1:6))
                    constr = [constr, eobs(k) == yobs(k) - y(k)]; % another way to enforce ydes?
                    
                end
                
                constr = [constr, u <= 100*ones(1,timeahead), u >= zeros(1,timeahead)];
                constr = [constr, y(end) == ydes];
                
                % cost function
                cost = norm(y - ydes*ones(1,timeahead+1));
                
                % call solvers
                optimize(constr,cost);
                
                % report results
                solu = value(u)
                soly = value(y);
                solx = value(x);
                sole = value(eobs);
                solc = value(cost)
                
%                 OPT.(pem{pm}).(tset{ts}).(order{od}).(timehorizon).solu = value(u);
%                 OPT.(pem{pm}).(tset{ts}).(order{od}).(timehorizon).soly = value(y);
%                 OPT.(pem{pm}).(tset{ts}).(order{od}).(timehorizon).solx = value(x);
%                 OPT.(pem{pm}).(tset{ts}).(order{od}).(timehorizon).sole = value(eobs);
%                 OPT.(pem{pm}).(tset{ts}).(order{od}).(timehorizon).solc = value(cost);
%                 
%                 cd '/Users/Aldis/Documents/MATLAB/ES100/Dataset4/Preprocessed/021319/Controltrial_031519_atrain'
%                 figure;
%                 % x axis build
%                 plotx = zeros(1,timeahead*2 + 2); % +2 for starting and end points
%                 ploty = zeros(1,timeahead*2 + 2);
%                 addedx = zeros(1);
%                 addedy = zeros(1,2);
%                 for r = 1:timeahead
%                     addx = repmat(r,1,2);
%                     addy = repmat(solu(r),1,2);
%                     addedx = horzcat(addedx, addx);
%                     addedy = horzcat(addedy, addy);
%                 end
%                 addx = timeahead+1;
%                 plotx = horzcat(addedx,addx);
%                 ploty = addedy;
%                 plot(plotx', ploty'); hold on; %... [0,0, solu(1),solu(1),solu(2),solu(2),solu(3),solu(3),... solu(4),solu(4),solu(5),solu(5)]
%                 plot(0:timeahead+1, [yobs(1), soly])
%                 xlabel('Minutes')
%                 ylabel('Temp ºC and % Open')
%                 legend('PositR (% open)','MPC-predicted roomC','Location','Northwest')
%                 ax = gca;
%                 ax.FontSize = 12;
%                 filesave = strcat(['controlsim_031519_',pem{pm},'_',tset{ts},'_',order{od},'_',num2str(timeahead),'min_trial',num2str(trialnum),'.png'])
%                 saveas(gcf,filesave)
%                 
%                 % with state space model alone (but using this control feedback...)
%                 % have a no-control (NC) and a with control state-space estimation
%                 states = zeros(ord,timeahead+1);
%                 states_nc = zeros(ord,timeahead+1);
%                 states(:,1) = xobs0;
%                 states_nc(:,1) = xobs0;
%                 roomtemp = zeros(1,timeahead+1);
%                 roomtemp_nc = zeros(1,timeahead+1);
%                 roomtemp(1) = C*states(:,1);
%                 roomtemp_nc(1) = C*states_nc(:,1);
%                 eobs = zeros(1,timeahead+1);
%                 eobs(1) = yobs(1)-roomtemp(1);
%                 for k = 2:timeahead+1
%                     % with control
%                     states(:,k) = A*states(:,k-1) + Bctrl*solu(k-1) + F*dobs0' + K*sole(k-1);
%                     % without control
% %                     states_nc(:,k) = A*states_nc(:,k) + B*[0,dobs0]' + K*eobs(k-1);
%                     
%                     roomtemp(k) = C*states(:,k) + sole(k);
% %                     roomtemp_nc(k) = C*states_nc(:,k) + eobs(k); % SHOULD BE STATES_NC!!!
% %                     eobs(k) = roomtemp_nc(k) - yobs(k); % but is this y the room temp changing due to control or no control?
%                 end
%                 
%                 % NEED CLARITY ON WHAT I'M TRYING TO SHOW USING STATE-SPACE NO CONTROL
%                 % BASELINE.. BECAUSE IT'S NOT THE LAB ITSELF...
%                 
%                 % CAN HAVE A FIGURE THAT SHOWS MODEL VS> ROOM BUT NOT MODEL VS MPC (BC
%                 % MPC USES MODEL)
%                 
%                 figure;
%                 plot(plotx', ploty'); hold on;
%                 plot(0:timeahead+1, [yobs(1), soly])
%                 plot(0:timeahead+1, [yobs(1),roomtemp],'o');
%                 plot(0:timeahead+1, [yobs(1),roomtemp_nc],'o');
%                 xlabel('Minutes')
%                 ylabel('Temp ºC and % Open')
%                 legend('PositR (% open)','MPC-predicted roomC','Model simulated control','Model simulated NO control','Location','Northwest')
%                 ax = gca;
%                 ax.FontSize = 12;
%                 filesave = strcat(['controlsim_031519_',pem{pm},'_',tset{ts},'_',order{od},'_',num2str(timeahead),'min_modsim_trial',num2str(trialnum),'.png'])
%                 saveas(gcf,filesave)
                
            end
%             close all;
        end
    end
end

solu
    
    % now we know it's a model problem
        % wasn't including K*e(t) where e(t) is observed-predicted roomC
    