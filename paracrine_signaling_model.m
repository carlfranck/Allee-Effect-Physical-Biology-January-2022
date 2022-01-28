%  SimBiology script for stochastic simulations of cell clustering
%  model by Carl Franck
%
%  Author: Igor Segota
%  Date:   03/03/2013

%  Note: Make sure the SimBiology toolbox is installed with MATLAB.
%  Note2: Get equations by getequations(Mobj)

%  Model
%  We will be simulating the following reactions
%  
%  Cell clustering:
%  1)   nn + nn -> nn + nc    (two single cells cluster with rate alpha)
%  2)   nn + nc -> nc + nc         (single cells and a cluster form another
%                              cluster)
%  3)   nc      -> nn         (clusters randomly desintegrate in single 
%                              cells)
%  Cell growth:
%  4)   nn      -> nn         (cell growth for single cells), rate gSlow
%  5)   nc      -> nc         (cell growth for cells in clusters)
%                                rate gFast
% 
%  Single cells grow slow, clustered cells grow fast
%
%  nn = number of cells not in a cluster
%  nc = number of cells in clusters (cluster sizes do not matter here;
%       as long as the cell is in 'some' cluster it is in a different
%       state than non-clustered cell nn)


% For paracrine model growth curve these numbers gave good fit with the
% mean growth curve for the bottle data (150RPM):
% alpha = (150/16)*1e-9;
% tau   = 2*1800; 
% tau is then 1 hour
% Tdslow = 18.3 hrs
% Tdfast = 9.3 hrs

% For vials let's try other numbers to get ~ 30 hour mean lag time
% (that we observe).


%% Build model


% Model parameters
alpha  = 0.45*1e-9;                 % clustering rate [1/s] default: 1e-9
tau    = 100*3600;             % cluster lifetime [s]
                               % measured 30 +- 20 min (0.5*3600)
                               %  ~ 1800 s
                               % measured on surface to be ~ ln(2)/2.5 hr-1
                               % which is about 20 minutes
n0     = 60;                   % initial number of single cells
                               %  default: 100
nc0    = 0;                    % initial number of cells in clusters
Tend   = 200*3600;             % Simulation end time in seconds

Tdslow = 18.3;                 % Doubling time in the slow (lag) phase [hr]
                               %  measured 18.3 hours
Tdfast = 11.3;                 % Doubling time in the fast (log) phase [hr]
                               %  measured 11.3 hours
gSlow  = log(2)/(Tdslow*3600); % Growth rate in slow phase [1/s]
gFast  = log(2)/(Tdfast*3600); % Growth rate in fast phase [1/s]

% Use this if we with to add stochastic growth rates:
% (normal dist. these are SDs, gSlow and gFast are then the means)
% gSlowSd= 0.01481979/3600;      % Slow phase growth rate StDev.
% gFastSd= 0.0208176/3600;       % Fast phase growth rate StDev.
gSlowSd= 0;
gFastSd= 0;
gDt    = 10;                   % Number of hours when resampling gamma

nRuns  = 100;                  % Run the model this many times to get
                               %  enough statistics.

% Model configuration
varRates  = 0;                 % Use variable rates
stopAtCx  = 0;                 % Stop the simulation first time the growth
                               % factor concentration reaches cx?
ccDecay   = 0;                 % Enable the growth factor random decay                               
ccStop    = 0;                 % Stop producing growth factor in clusters
                               %  when cc>cx


%% Solve the SSA model nRuns times

% I used this for histogram in the paper

% Model configuration
varRates  = 0;                 % Use variable rates
stopAtCx  = 1;                 % Stop the simulation first time the growth
                               % factor concentration reaches cx?
ccDecay   = 0;                 % Enable the growth factor random decay                               
ccStop    = 1;                 % Stop producing growth factor in clusters
                               %  when cc>cx

nRuns        = 50;
lagTimesHrsSSA = zeros(1, nRuns);

for i=1:nRuns 
    sbiomodel_clustering;
    
    cs = getconfigset(Mobj);
    set(cs, 'StopTime', Tend);
    set(cs, 'SolverType', 'ssa');
    cs.SolverOptions.LogDecimation = 1000;
    [t_ssa, x_ssa, names]   = sbiosimulate(Mobj);
    
    % Find the lag time
    iEnd = length(x_ssa(:,1));     % Index of the last point
    nEnd = x_ssa(iEnd,1)+x_ssa(iEnd,2);  % total density at last point
    tLag = Tend - (1/gFast)*log(nEnd/n0);
    tLagHrs = tLag/3600;

    lagTimesHrsSSA(i) = tLagHrs;

    % Progress
    disp(i/nRuns)
end

mean(lagTimesHrsSSA)
std(lagTimesHrsSSA)
hist(lagTimesHrsSSA)

% 
% FIG0 = figure;
% set(gcf, 'color', 'white');
% semilogy(t_ssa/3600, x_ssa(:,([1 2 3])));
% title('Paracrine signaling SSA model');
% xlabel('time (hours)');
% ylabel('amount per 1 ml (molecules)');
% legend(names([1 2 3]), 'Location', 'NorthEastOutside');
