% hypoxia settings:
% ischeaemia with 0.0001 scalar to RhoREN and RhoREF (Electron carrier
% conc) (0.01% oxygen)


%% new starting states for next simulation of 1 hz for 10 minutes?

paramIX.bcl = 1000; % basic cycle length in ms
paramIX.verbose = true; % printing numbers of beats simulated.
paramIX.flag_ode =  1;
paramIX.stimDur = 1; % no stimulation to get steady state
paramIX.nao = 137;
paramIX.ko = 5.4;
paramIX.cao = 1.8;
paramIX.model = @model_Torord_withmito_v6;

options = []; % parameters for ode15s - usually empty
beatsIX = 600; % number of beats
ignoreFirst = 0; % this many beats at the start of the simulations are ignored when extracting the structure of simulation outputs (i.e., beats - 1 keeps the last beat).

% make an array of parameter structures, 3 for one of each model
paramsIX(1:4) = paramIX;

% and then each is assigned a different model
paramsIX(1).model = @model_Torord_withmito_v6;
paramsIX(2).model = @model_Torord_withmito_v6_normoxia; 
paramsIX(3).model = @model_Torord_withmito_v6_modhypoxia; 
paramsIX(4).model = @model_Torord_withmito_v6_sevhypoxia; 

%% Simulation and extraction of outputs

parfor i=1:4
    X0_0{i} = getStartingState_withmito_v6('Torord_endo'); 
    [time_IX{i}, X_IX{i}] = modelRunner_withmito_v6(X0_0{i}, options, paramsIX(i), beatsIX, ignoreFirst);
    currentsIX{i} = getCurrentsStructure_withmito_v6(time_IX{i}, X_IX{i}, paramsIX(i), 0);
end


%% figure 0.1hz ATP

control_ATPi_cyto = mean([max(X_IX{1}{1,1}(:,56)), min(X_IX{1}{1,1}(:,56))]);
normox_ATPi_cyto = mean([max(X_IX{2}{1,1}(:,56)), min(X_IX{2}{1,1}(:,56))]);
modhypoxia_ATPi_cyto = mean([max(X_IX{3}{1,1}(:,56)), min(X_IX{3}{1,1}(:,56))]);
sevhypoxia_ATPi_cyto = mean([max(X_IX{4}{1,1}(:,56)), min(X_IX{4}{1,1}(:,56))]);


ATPi_cyto = [control_ATPi_cyto, normox_ATPi_cyto, modhypoxia_ATPi_cyto, sevhypoxia_ATPi_cyto];

control_CRP_cyto = mean([max(X_IX{1}{1,1}(:,59)), min(X_IX{1}{1,1}(:,59))]);
normox_CRP_cyto = mean([max(X_IX{2}{1,1}(:,59)), min(X_IX{2}{1,1}(:,59))]);
modhypoxia_CRP_cyto = mean([max(X_IX{3}{1,1}(:,59)), min(X_IX{3}{1,1}(:,59))]);
sevhypoxia_CRP_cyto = mean([max(X_IX{4}{1,1}(:,59)), min(X_IX{4}{1,1}(:,59))]);

CRP_cyto = [control_CRP_cyto, normox_CRP_cyto, modhypoxia_CRP_cyto, sevhypoxia_CRP_cyto];

vals = [ATPi_cyto ; CRP_cyto];
x = [1 2 3 4];

hypoxia = categorical(["Control", "Normoxia", "Moderate Hypoxia", "Severe Hypoxia"]);
hypoxia = reordercats(hypoxia, ["Control", "Normoxia", "Moderate Hypoxia", "Severe Hypoxia"]);

figure(1);
b = bar(hypoxia, vals);
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [1 1 1];
ylabel("ATP or CrP cytoplasm (mM)")

%% change time to s
for i=1:4
    for k = 1:beatsIX
        time_IX_s{1,i}{k, 1} = time_IX{1,i}{k, 1} ./ 1000;
        time_IX_s{1,i}{k, 1} = time_IX_s{1,i}{k, 1} + 1*(k-1);
    end
end

%% plot severe hypoxia as ischaemia equivalent
figure(2);
for k = 1:beatsIX
    plot(time_IX_s{1,4}{k, 1}, X_IX{4}{k,1}(:,57), "black")
    plot(time_IX_s{1,4}{k, 1}, X_IX{4}{k,1}(:,59), "cyan")
    plot(time_IX_s{1,4}{k, 1}, X_IX{4}{k,1}(:,2), "red")
    hold on
end
xlabel("time (s)")
ylabel("Substrate concentration (mM)")

qw{1} = plot(nan, 'k-');
qw{2} = plot(nan, 'c-');
qw{3} = plot(nan, 'r-');
legend([qw{:}], {'ATPi_{cyto}','CrPi_{cyto}','Nai'}, 'location', 'best')
