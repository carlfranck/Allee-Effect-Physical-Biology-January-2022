%  SimBiology script for stochastic simulations of cell clustering
%  model; endocrine signaling.
%
%  Author: Igor Segota
%  Date:   02/07/2013

%  Note: Make sure the SimBiology toolbox is installed with MATLAB.

%  Model
%  We will be simulating the following reactions
%  
%  Cell clustering:
%  1)   nn + nn -> nc         (two single cells cluster with rate alpha)
%  2)   nn + nc -> nc         (single cells and a cluster form another
%                              cluster)
%  3)   nc      -> nn         (clusters randomly desintegrate in single 
%                              cells)
%  Cell growth:
%  4)   nn      -> nn         (cell growth for single cells)
%  5)   nc      -> nc         (cell growth for cells in clusters)
% 
%  Growth factor dynamics:
%  6)   nc      -> nc + cc
%  nn = number of cells not in a cluster
%  nc = number of cells in clusters (cluster sizes do not matter here;
%       as long as the cell is in 'some' cluster it is in a different
%       state than non-clustered cell nn)

% Model options
clu    = 1;                    % enable(1) or disable(1) clustering
grow   = 1;                    % enable(1) or disable(1) cell growth

% Model parameters
alpha  = (2e-9);               % clustering rate [1/s] default: 1e-9
tau    = 0.5*3600;             % cluster lifetime [s]; 
                               % Measured to be around 30 minutes based on
                               % surface growth data.
                               % See: "times between dissociation events into two things.txt"
n0     = 60;             % initial number of single cells
                               %  default: 100
nc0    = 0;                    % initial number of cells in clusters
c0     = 0;                    % initial growth factor number
Tend   = 200*3600;             % Simulation end time in seconds
                               %  default: 200 hr *3600 s/hr
cx     = 7.5*10^(9);             % Growth factor number at transition
                               %  default: 7e8
                               % Calculation for Kd = 0.5 nM is performed
                               % in the file 
                               % Growth factor concentration at the transition for cluster-based endocrine and paracrine models.txt
                               % and gives us ~ 7.5 * 10^12
Tdslow = 19.5;                 % Doubling time in the slow (lag) phase [hr]
Tdfast = 11.3;                 % Doubling time in the fast (log) phase [hr]
                               % Tdfast = 9.5 hrs for vialrs, 11.3 hrs for
                               % bottles
gSlow  = log(2)/(Tdslow*3600); % Growth rate in slow phase [1/s]
gFast  = log(2)/(Tdfast*3600); % Growth rate in fast phase [1/s]
nu     = 2000;                 % Growth factor secretion rate per cell
                               %  [molecules/(cell*s)] default: 9000

% Dimensionless model parameters   
%  I'm not sure if I want to use these.
TTend  = alpha*Tend;
alpTau = alpha*tau;
gpSlow = gSlow/alpha;
gpFast = gFast/alpha;
alphap = 1;

% This is where we create the main model
Mobj    = sbiomodel('cellClustering');
compObj = addcompartment(Mobj, 'contents');

% Clustering reactions
% Robj are reactions within the main model
% Kobj are kinetic laws, I usually go with MassAction kinetics
% Pobj are parameters


if (clu==1)
    Robj1   = addreaction(Mobj, '2 nn -> nc + nn');
    Kobj1   = addkineticlaw(Robj1, 'MassAction');
    Pobj1   = addparameter(Mobj, 'alpha', alpha);
    set(Kobj1, 'ParameterVariableNames', 'alpha');

    Robj2   = addreaction(Mobj, 'nn + nc -> 2 nc');
    Kobj2   = addkineticlaw(Robj2, 'MassAction');
    set(Kobj2, 'ParameterVariableNames', 'alpha');

    Robj3   = addreaction(Mobj, 'nc -> nn');
    Kobj3   = addkineticlaw(Robj3, 'MassAction');
    Pobj3   = addparameter(Kobj3, 'inverseTau', 1/tau);
    set(Kobj3, 'ParameterVariableNames', 'inverseTau');
end 

% Growth reactions
if (grow==1)
    Robj4   = addreaction(Mobj, 'nn -> 2 nn');
    Kobj4   = addkineticlaw(Robj4, 'MassAction');
    Robj5   = addreaction(Mobj, 'nc -> 2 nc');
    Kobj5   = addkineticlaw(Robj5, 'MassAction');
    Pobj45  = addparameter(Mobj, 'gamma', gSlow, 'ConstantValue', false);
    
    set(Kobj4, 'ParameterVariableNames', 'gamma');
    set(Kobj5, 'ParameterVariableNames', 'gamma');

    % Use this event for deterministic simulation to get the entire
    % growth curve. 
    Mobj.addevent(['cc >= ' num2str(cx)], ...
        {['gamma = ' num2str(gFast)]});
    
    % Use this event for stochastic simulation, since we do not need
    % to keep tracking growth once we're in fast growth phase.
    % This sets all chemical species to 0 to speed up the calculation.
    %Mobj.addevent(['cc >= ' num2str(cx)], ...
    %    {'nn = 0', 'nc = 0', 'cc = 0', 'nu = 0'});
    
    % Put gFast instead of 0 here if we want to see the log phase dynamics
    % Stop producing the growth factor in fast growing (log) phase
end

% Growth factor production reaction
Robj6   = addreaction(Mobj, 'nc -> nc + cc');
Kobj6   = addkineticlaw(Robj6, 'MassAction');
Pobj6   = addparameter(Mobj, 'nu', nu, 'ConstantValue', false);
set(Kobj6, 'ParameterVariableNames', 'nu');


% Set initial conditions
Sobj1 = sbioselect(Mobj, 'Type', 'species', 'Name', 'nn');
set(Sobj1, 'InitialAmount', n0);
Sobj2 = sbioselect(Mobj, 'Type', 'species', 'Name', 'nc');
set(Sobj2, 'InitialAmount', nc0);
Sobj3 = sbioselect(Mobj, 'Type', 'species', 'Name', 'cc');
set(Sobj3, 'InitialAmount', c0);

% Solve the ODE (mean field) model
cs = getconfigset(Mobj);
set(cs, 'StopTime', Tend);
%set(cs, 'SolverType', 'ssa');
%cs.SolverOptions.LogDecimation = 100000;
%cs.SolverOptions.AbsoluteTolerance = 1e-30;
% cs.SolverOptions.RelativeTolerance = 1e-12;
%cs.SolverOptions.MaxStep = 1;
[t_ode, x_ode, names] = sbiosimulate(Mobj);















FIG0 = figure;
set(gcf, 'color', 'white');
semilogy(t_ode/3600, (x_ode(:,(1))+x_ode(:,(2)))/25 );
title('Cell clustering ODE model');
xlabel('time (hours)');
ylabel('cell density (cells/ml)');
legend(names, 'Location', 'NorthEastOutside');
%axis([0, Tend, 1e-2, 1e6])

semilogy(t_ode/3600, (x_ode(:,(3)))/25 );

%plot(t_ode, g_ode)

% How to check what solver is currently selected?
%  getconfigset(Mobj)
% How to change the solver?
%  cs = getconfigset(Mobj)
%  set(cs, 'SolverType', 'ssa')
%  cs
% Change the value of log decimation
%  cs.SolverOptions.LogDecimation = 10;
%  Increasing this setting lets us record fewer data points and decrease
%  run time.
% [t_ssa, x_ssa] = sbiosimulate(Mobj, cs)

% dlmwrite()
