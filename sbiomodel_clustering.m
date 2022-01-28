                              
Mobj    = sbiomodel('cellClustering');
compObj = addcompartment(Mobj, 'contents');

% Growth reactions
Robj4   = addreaction(Mobj, 'nn -> 2 nn');
Kobj4   = addkineticlaw(Robj4, 'MassAction');
Robj5   = addreaction(Mobj, 'nc -> 2 nc');
Kobj5   = addkineticlaw(Robj5, 'MassAction');
Pobj4   = addparameter(Mobj, 'gammaSlow', gSlow);
Pobj5   = addparameter(Mobj, 'gammaFast', gFast);
PobjF   = addparameter(Mobj, 'gammaBool', 0, 'ConstantValue', false);
set(Kobj4, 'ParameterVariableNames', {'gammaSlow'});
set(Kobj5, 'ParameterVariableNames', {'gammaFast'});

% Clustering reactions
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

% Growth factor production reaction
%Robj6   = addreaction(Mobj, 'nc -> nc + cc');
%Kobj6   = addkineticlaw(Robj6, 'MassAction');
%Pobj6   = addparameter(Mobj, 'nu', nu, 'ConstantValue', false);
%set(Kobj6, 'ParameterVariableNames', 'nu');
%Robj6a  = addreaction(Mobj, 'cc -> null');
%Kobj6a  = addkineticlaw(Robj6a, 'MassAction');
%Pobj6a  = addparameter(Kobj6a, 'inverseTauC', 1/tauC);
%set(Kobj6a, 'ParameterVariableNames', 'inverseTauC');


% If we allow for the rate parameters to change in time then

if (varRates == 1)
    rateChangeDt = gDt*3600;

    for tt = 0:rateChangeDt:Tend
        Mobj.addevent(['time >= ' num2str(tt) ' && cc < ' num2str(cx)], ...
            ['gamma = ' num2str(gSlow + gSlowSd*randn(1))]);
        Mobj.addevent(['time >= ' num2str(tt) ' && cc >= ' num2str(cx)], ...
            ['gamma = ' num2str(gFast + gFastSd*randn(1))]);
    end
end


%
% No growth factor secretion here so this is commented out.
%


% set gammaBool to 1 when we cross the cx first time. Later we'll use
% this for stopping cell growth
%Mobj.addevent(['cc >= ' num2str(cx) ' && gammaBool == 0'], ...
%    'gammaBool = 1');

%if (stopAtCx == 1)
%    Mobj.addevent(['cc >= ' num2str(cx)], ...
%        {'nn = 0', 'nu = 0', 'nc = 0', 'gamma = 0', 'cc = 0'});
%else
%    if (ccStop == 1)
%        Mobj.addevent(['cc >= ' num2str(cx)], ...
%            {['gamma = ' num2str(gFast)], ['nu = ' num2str(nu/1e5)]});
%        Mobj.addevent(['cc <= ' num2str(cx/10)], ...
%            {['gamma = ' num2str(gSlow)], ['nu = ' num2str(nu)]});
%        %Mobj.addevent(['cc <= ' num2str(cx/10)], ...
%        %    'nu = 0.09');
%        % ['nu = ' num2str(nu)]
%    else
%        Mobj.addevent(['cc >= ' num2str(cx)], ...
%            ['gamma = ' num2str(gFast)]);
%        Mobj.addevent(['cc < ' num2str(cx2)], ...
%            ['gamma = ' num2str(gSlow)]);
 %   end
%end


% Set initial conditions
Sobj1 = sbioselect(Mobj, 'Type', 'species', 'Name', 'nn');
set(Sobj1, 'InitialAmount', n0);
Sobj2 = sbioselect(Mobj, 'Type', 'species', 'Name', 'nc');
set(Sobj2, 'InitialAmount', nc0);
