% Simulation to recover fig3 on Cortassa 2006
% Simulation protocol is as follows:
% no stim for 800s
% 0.25 hz for 800s
% 0.5 hz for 200s
% 1.0 hz for 200s
% 1.5 hz for 200s
% 2.0 hz for 200s

% hypoxia scalar to rhoREN and rhoREF:
% normoxia: 0.21
% moderate hypoxia: 0.005
% severe hypoxia: 0.001


% no stim for 800s
paramI.bcl = 4000; % basic cycle length in ms
paramI.verbose = true; % printing numbers of beats simulated.
paramI.flag_ode =  1;
paramI.stimDur = 0; 
paramI.nao = 137;
paramI.ko = 5.4;
paramI.cao = 1.8;
paramI.model = @model_Torord_withmito_v6_VnaCa_unscaled;

options = []; % parameters for ode15s - usually empty
beatsI = 200; % number of beats
ignoreFirst = 0; % this many beats at the start of the simulations are ignored when extracting the structure of simulation outputs (i.e., beats - 1 keeps the last beat).

% make an array of parameter structures, 3 for one of each model
paramsI(1:4) = paramI;

% and then each is assigned a different model
paramsI(1).model = @model_Torord_withmito_v6_VnaCa_unscaled;
paramsI(2).model = @model_Torord_withmito_v6_normoxia; 
paramsI(3).model = @model_Torord_withmito_v6_modhypoxia; 
paramsI(4).model = @model_Torord_withmito_v6_sevhypoxia; 

%% Simulation and extraction of outputs

parfor i=1:4
    X0_I = getStartingState_withmito_v6('Torord_endo'); % starting state - can be also Torord_mid or Torord_epi for midmyocardial or epicardial cells respectively.
    [time_I{i}, X_I{i}] = modelRunner_withmito_v6(X0_I, options, paramsI(i), beatsI, ignoreFirst);
    currentsI{i} = getCurrentsStructure_withmito_v6(time_I{i}, X_I{i}, paramsI(i), 0);
end

%% new starting states for next simulation of 0.25 hz for 800s

paramII.bcl = 4000; % basic cycle length in ms
paramII.verbose = true; % printing numbers of beats simulated.
paramII.flag_ode =  1;
paramII.stimDur = 1; % no stimulation to get steady state
paramII.nao = 137;
paramII.ko = 5.4;
paramII.cao = 1.8;
paramII.model = @model_Torord_withmito_v6_VnaCa_unscaled;

options = []; % parameters for ode15s - usually empty
beatsII = 200; % number of beats
ignoreFirst = 0; % this many beats at the start of the simulations are ignored when extracting the structure of simulation outputs (i.e., beats - 1 keeps the last beat).

% make an array of parameter structures, 3 for one of each model
paramsII(1:4) = paramII;

% and then each is assigned a different model
paramsII(1).model = @model_Torord_withmito_v6_VnaCa_unscaled;
paramsII(2).model = @model_Torord_withmito_v6_normoxia; 
paramsII(3).model = @model_Torord_withmito_v6_modhypoxia; 
paramsII(4).model = @model_Torord_withmito_v6_sevhypoxia; 

%% Simulation and extraction of outputs

parfor i=1:4
    X0_II{i} = X_I{i}{end,1}(end,:);
    [time_II{i}, X_II{i}] = modelRunner_withmito_v6(X0_II{i}, options, paramsII(i), beatsII, ignoreFirst);
    currentsII{i} = getCurrentsStructure_withmito_v6(time_II{i}, X_II{i}, paramsII(i), 0);
end


%% new starting states for next simulation of 1 hz for 10 minutes?

paramIX.bcl = 1000; % basic cycle length in ms
paramIX.verbose = true; % printing numbers of beats simulated.
paramIX.flag_ode =  1;
paramIX.stimDur = 1; % no stimulation to get steady state
paramIX.nao = 137;
paramIX.ko = 5.4;
paramIX.cao = 1.8;
paramIX.model = @model_Torord_withmito_v6_VnaCa_unscaled;

options = []; % parameters for ode15s - usually empty
beatsIX = 600; % number of beats
ignoreFirst = beatsIX - 1; % this many beats at the start of the simulations are ignored when extracting the structure of simulation outputs (i.e., beats - 1 keeps the last beat).

% make an array of parameter structures, 3 for one of each model
paramsIX(1:4) = paramIX;

% and then each is assigned a different model
paramsIX(1).model = @model_Torord_withmito_v6_VnaCa_unscaled;
paramsIX(2).model = @model_Torord_withmito_v6_normoxia; 
paramsIX(3).model = @model_Torord_withmito_v6_modhypoxia; 
paramsIX(4).model = @model_Torord_withmito_v6_sevhypoxia; 

%% Simulation and extraction of outputs

parfor i=1:4
    X0_0{i} = X_II{1}{end,1}(end,:);
    [time_IX{i}, X_IX{i}] = modelRunner_withmito_v6(X0_0{i}, options, paramsIX(i), beatsIX, ignoreFirst);
    currentsIX{i} = getCurrentsStructure_withmito_v6(time_IX{i}, X_IX{i}, paramsIX(i), 0);
end


%% figure 1hz ATP

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

%% ATP CrP normalised

ATPi_cyto_norm = [normox_ATPi_cyto modhypoxia_ATPi_cyto sevhypoxia_ATPi_cyto]./ normox_ATPi_cyto * 100;
CRP_cyto_norm = [normox_CRP_cyto modhypoxia_CRP_cyto sevhypoxia_CRP_cyto]./ normox_CRP_cyto * 100;

data_ATPi_norm = [100 82.7 86.5];
data_CRP_norm = [100 111 63.3];

sims_data_norm = [data_ATPi_norm; ATPi_cyto_norm; data_CRP_norm; CRP_cyto_norm];

hypoxia_norm = categorical([ "Normoxia", "Moderate Hypoxia", "Severe Hypoxia"]);
hypoxia_norm = reordercats(hypoxia_norm, ["Normoxia", "Moderate Hypoxia", "Severe Hypoxia"]);

figure(2);
b = bar(hypoxia_norm, sims_data_norm);
set(b(1),"FaceColor", [0 0.4470 0.7410])
set(b(2),"FaceColor", [1 1 1],'EdgeColor', [0 0.4470 0.7410], 'LineWidth', 2)
set(b(3),"FaceColor", [0.9290 0.6940 0.1250])
set(b(4),"FaceColor", [1 1 1],'EdgeColor', [0.9290 0.6940 0.1250], 'LineWidth', 2)
ylabel("ATP or CrP Normalised (%)");
legend([b(1) b(2) b(3) b(4)], "Measured ATP", "Simulated ATP", "Measured CrP", "Simulated CrP")
set(gcf,'color','w');

%% figure for oxygen consumption rate

control_VNO = mean([max(currentsIX{1,1}.VNO), min(currentsIX{1,1}.VNO)])*1000*60;
normox_VNO = mean([max(currentsIX{1,2}.VNO), min(currentsIX{1,2}.VNO)])*1000*60;
midhypox_VNO = mean([max(currentsIX{1,3}.VNO), min(currentsIX{1,3}.VNO)])*1000*60;
sevhypox_VNO = mean([max(currentsIX{1,4}.VNO), min(currentsIX{1,4}.VNO)])*1000*60;

VNO_plot = [control_VNO, normox_VNO, midhypox_VNO, sevhypox_VNO];

figure(3);
hypoxia = categorical(["Control", "Normoxia", "Moderate Hypoxia", "Severe Hypoxia"]);
hypoxia = reordercats(hypoxia, ["Control", "Normoxia", "Moderate Hypoxia", "Severe Hypoxia"]);
bar(hypoxia, VNO_plot, "black")
ylabel("Oxygen Consumption Rate (nM.min^{-1})")

%% figure for oxygen consumption rate normalised

VNO_plot_normalised = [normox_VNO, midhypox_VNO, sevhypox_VNO] / normox_VNO * 100;

VNO_plot_data = [100, 78.47, 0];

VNO_plot_all = [VNO_plot_data; VNO_plot_normalised];

hypoxia_norm = categorical([ "Normoxia", "Moderate Hypoxia", "Severe Hypoxia"]);
hypoxia_norm = reordercats(hypoxia_norm, ["Normoxia", "Moderate Hypoxia", "Severe Hypoxia"]);

figure(4);
b = bar(hypoxia_norm, VNO_plot_all);
set(b(1),"FaceColor", [0 0 0])
set(b(2),"FaceColor", [1 1 1],'EdgeColor', [0 0 0], 'LineWidth', 2)
legend([b(1) b(2)], "Measured", "Simulated")
ylabel("Oxygen Consumption Rate Normalised (%)")
ylim([0 110])
set(gcf,'color','w');
