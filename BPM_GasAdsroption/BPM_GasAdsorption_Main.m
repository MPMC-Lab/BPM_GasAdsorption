%% General Code Description of BPM of gas adsorption and Author Introduction
%======================================================================================================================
%> @ Code Description:
%> @ File        PI_GasAdsorption.m
%> @ Brief       Advanced parameter identification module based on Bayesian inference for gas adsorption processes.
%> @ Details     PI_GasAdsorption serves as the main program for parameter identification in gas adsorption models,
%>               encompassing a holistic approach that integrates experimental data processing with Bayesian inference.
%>
%>               The process is characterized by the following steps:
%>               (1) Parameter Establishment: Initiation by defining parameters based on physical attributes of gases and adsorbents.
%>               (2) Bayesian Inference Optimization: Alignment of Bayesian inference settings with material properties and experimental conditions,
%>                   including optimization of steps, samples, and calibration of model predictions against actual measurements.
%>               (3) Iterative Refinement and Analysis: Ensures efficient and effective Bayesian inference iterations, leading to a comprehensive
%>                   analysis where estimated parameters from the model are juxtaposed with actual experimental data, often visualized through graphs.
%>               (4) Module for Bayesian Parameter Identification: Facilitates the entire process of parameter identification,
%>                   from setting initial values and options for Bayesian inference to numerically solving the model with estimated parameters
%>                   and optimizing these parameters through a rigorous Bayesian analysis.
%>
%> @ Author Introduction
%>              - Yesol Hyun (yesol2@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University
%>              - Geunwoo Oh (gwoh@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University
%>              - Jung-Il Choi (jic@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University
%>
%> @ Date        November 2023
%> @ Version     1.0
%> @ Copyright   Copyright (c) 2023-2023 Yesol Hyun, Geunwoo Oh, and Jung-Il Choi, Yonsei University, All rights reserved.
%> @ License     This project is released under the terms of the MIT License (see LICENSE file).
%======================================================================================================================
close all;clear all; clc; warning off;
rng(100,'twister')
uqlab
%% Experimental Data Preparation and Initial Analysis for Parameter Identification
%======================================================================================================================
%> @details     This section of the PI_GasAdsorption.m module encompasses the preparation of experimental data and the initial setup for Bayesian analysis in gas adsorption models. The process includes:

%>              (1) Data Loading and Processing:
%>                  - Loading and calculating properties such as flow rates, bed length, particle diameter, bed volume, bulk density, and superficial velocity.
%>                  - Reading and processing experimental data from CSV files, including time and concentration data.

%>              (2) Data Analysis and Preparation for Bayesian Inference:
%>                  - Performing linear interpolation, numerical integration, and normalization of data.
%>                  - Preparing data for Bayesian inference by interpolating and smoothing concentration data.
%>                  - Compiling differential concentration rate data and normalized time-to-breakthrough data.

%>              (3) Bayesian Analysis Setup:
%>                  - Defining the forward model for the Bayesian analysis using a custom function.
%>                  - Setting prior distributions for various parameters based on expected ranges.
%>                  - Configuring the Bayesian inversion solver, including the selection of a sampler and defining the number of chains and steps.
%>                  - Preparing the discrepancy model and options for the Bayesian analysis.

%>              This comprehensive approach ensures that all the necessary steps are taken to accurately prepare the experimental data and set up the Bayesian analysis for effective parameter identification in gas adsorption models.
%======================================================================================================================
% Load measured data
Q_read = [21.1, 42.3, 42.3, 42.3, 63.5]; % Flow rates [sccm]
L_read = [0.01, 0.01, 0.02, 0.03, 0.03] * 100; % Bed length [cm]
dz = 0.01 / 50; % Grid size for the spatial discretization [m]
NumCurve = length(L_read); % Number of experimental data curves
L = L_read / 100; % Convert bed length from cm to m
cFeed = 0.002; % Feed concentration [mol/m^3]
ngrid = L / dz + 1; % Number of grid points [m]
mass = [79.8, 79.8, 157.1, 241.5, 241.3] * 1e-6; % Mass of adsorbent [kg]
Q = Q_read / 60 * (1e-2)^3; % Convert flow rates from cm^3/min to m^3/s
Diam = 0.004; % Diameter of the column [m]
A = (Diam / 2)^2 * pi; % Cross-sectional area of the column [m^2]
V_ads = A .* L; % Volume of the adsorbent [m^3]
rho_b = mass ./ V_ads; % Bulk density of the adsorbent [kg/m^3]
v = Q / A; % Superficial velocity [m/s]
dp = 0.00023; % Particle diameter [m]
mu = 1.8e-5; % Dynamic viscosity of the gas [Pa.s]
rho_g = 0.05795; % Density of the gas [kg/m^3]
epsilon = 0.39; % Bed porosity
rho = rho_b / epsilon; % Adjusted density [kg/m^3]
dt = 20; % Time step [s]
nz = ngrid; % Alias for number of grid points
us = v; % Superficial velocity alias [m/s]
Ndata = 30; % Number of data points for Bayesian analysis

ModelNum =1;
if ModelNum ==1
    ModelName = 'Langmuir';
elseif ModelNum==2
    filename = 'Linear';
elseif  ModelNum==3
    filename = 'Toth';
end

% Processing each experimental data curve
Max_C= 0.95; 
for i = 1:NumCurve
    % Load experimental data from CSV files
    filename = sprintf('ASZM-TEDA+GB(%.1fsccm,%dcm).csv', Q_read(i), L_read(i));
    exp_inf = csvread(filename, 1, 0);
    
    % Extract and preprocess time and concentration data
    t_exp_tmp = exp_inf(:, 1);
    idx_tmax = find(t_exp_tmp == max(t_exp_tmp), 1, 'first');
    t_data_tmp = t_exp_tmp(1:idx_tmax) * 60; % Convert time to seconds
    c_tmp = exp_inf(1:idx_tmax, 2);
    
    % Interpolate concentration data for uniform time distribution
    t_data = linspace(t_data_tmp(1), max(t_data_tmp), 20 * Ndata);
    c_data = pchip(t_data_tmp, c_tmp, t_data);
    
    % Calculate integral results for Bayesian inference
    max_c_limit = 0.95 * max(c_data);
    interp_func = @(t) interp1(t_data, 1 - c_data, t, 'linear', 'extrap');
    integration_limits = [0, max(t_data(1:end-1)), max_c_limit];
    integral_result(i) = cFeed * Q(i) * integral(interp_func, integration_limits(1), integration_limits(2)) / mass(i);
    
    % Further data processing for Bayesian inference
    idx1 = find(c_data <= 0.1, 1, 'last');
    idx2 = find(c_data > 0 & c_data <= Max_C, 1, 'last');
    
    tf(i)= max(t_data);
    
    t_interp = linspace(t_data(idx1), t_data(idx2), Ndata);
    c_interp = pchip(t_data(idx1:idx2), c_data(idx1:idx2) * cFeed, t_interp);
    
    % Numerical differentiation for rate of change in concentration
    dt_tmp = t_interp(2) - t_interp(1);
    dc_dt_new = diff(smooth(t_interp, c_interp, 0.9, 'lowess')) / dt_tmp;
    max_dcdt(i) = max(dc_dt_new);
    data_dcdt = dc_dt_new' / max_dcdt(i);
    
    % Determining specific times for a given concentration level
    specific_c = 0.01;
    specific_t = pchip(c_data(idx1:length(c_data)), t_data(idx1:length(c_data)) / 60, specific_c);
    data_tbr_tmp(i) = specific_t;
    
    % Compiling all processed data for Bayesian analysis
    if i == 1
        dcdt_tmp = data_dcdt;
        t_exp = t_interp;
        c_exp = c_interp;
    else
        dcdt_tmp = [dcdt_tmp, data_dcdt];
        t_exp = [t_exp, t_interp];
        c_exp = [c_exp, c_interp];
    end
end

% Finalizing the output data for Bayesian analysis : Normalize data
max_tbr = max(data_tbr_tmp);
data_tbr = data_tbr_tmp / max(data_tbr_tmp);
myData.y = [((Ndata - 1) / Ndata) * dcdt_tmp, (1 / Ndata) * data_tbr];
myData.Name = 'Exp data';

% Setting up the forward model for Bayesian analysis
ModelOpts.mHandle = @(par) uq_GasAdsorption(par, t_exp, c_exp, dt, ngrid, dz, cFeed, epsilon, tf, v,ModelNum,rho, rho_b, rho_g, dp, mu, NumCurve, max_dcdt, max_tbr, Ndata,Max_C);
ModelOpts.isVectorized = true;
myForwardModel = uq_createModel(ModelOpts);
BayesOpts.ForwardModel = myForwardModel;

% Defining prior distributions for Bayesian parameters
PriorOpts.Marginals(1).Name = 'b1';
PriorOpts.Marginals(1).Type = 'Uniform';
PriorOpts.Marginals(1).Parameters = [1 5];

PriorOpts.Marginals(2).Name = 'qm';
PriorOpts.Marginals(2).Type = 'Uniform';
PriorOpts.Marginals(2).Parameters = [-2 0];

PriorOpts.Marginals(3).Name = 'Dm';
PriorOpts.Marginals(3).Type = 'Uniform';
PriorOpts.Marginals(3).Parameters = [-7 -5];

PriorOpts.Marginals(4).Name = 'De';
PriorOpts.Marginals(4).Type = 'Uniform';
PriorOpts.Marginals(4).Parameters = [-7 -5];

% Creating the prior distribution as an input object
myPriorDist = uq_createInput(PriorOpts);
BayesOpts.Prior = myPriorDist;

% Setting up the discrepancy model for the Bayesian analysis
SigmaOpts.Marginals.Name = 'sigma2';
SigmaOpts.Marginals.Type = 'Uniform';
SigmaOpts.Marginals.Parameters = [0 (0.05)^2];
mySigmaDist = uq_createInput(SigmaOpts);
DiscrepancyOptsUnknownDisc.Type = 'Gaussian';
DiscrepancyOptsUnknownDisc.Prior = mySigmaDist;
BayesOpts.Discrepancy = DiscrepancyOptsUnknownDisc;

% Configuring the Bayesian inversion solver and sampler
Solver.Type = 'MCMC';
Solver.MCMC.Sampler = 'AIES';
Solver.MCMC.NChains = 100;
Solver.MCMC.Visualize.Interval = 10;
Solver.MCMC.Steps =1000;
Solver.MCMC.Visualize.Parameters = [1 2 3 4];

% Finalizing Bayesian inversion setup
BayesOpts.Type = 'Inversion';
BayesOpts.Data = myData;
BayesOpts.Solver = Solver;
%% Execution of Bayesian Analysis and Forward Modeling Using IUQ Results
%======================================================================================================================
%> @details     This section performs the Bayesian analysis for gas adsorption parameter identification and conducts forward modeling using the inferred parameters.
%>              The process includes:

%>              (1) Bayesian Analysis Execution:
%>                  - Perform Bayesian analysis using the specified options.
%>
%>              (2) Post-Processing of Bayesian Results:
%>                  - Filter out MCMC chains with low acceptance rates.
%>                  - Display and analyze the acceptance rates of the chains.
%>
%>              (3) Results Reporting and Visualization:
%>                  - Print a detailed report of the Bayesian analysis results.
%>                  - Display results including prior and posterior distributions.
%>
%>              (4) Forward Modeling Based on Inferred Parameters:
%>                  - Utilize inferred parameters to run forward models.
%>                  - Compare modeled data with experimental data and visualize the results.
%======================================================================================================================
% Execution of Bayesian Analysis
tic
myBayesianAnalysis_Step = uq_createAnalysis(BayesOpts); % Create Bayesian analysis
uq_time = toc; % Record the time taken for the analysis

% Filtering MCMC chains with low acceptance rates
uq_display(myBayesianAnalysis_Step, 'acceptance', true); % Display acceptance rates
acceptance = myBayesianAnalysis_Step.Results.Acceptance; % Extract acceptance rates
[~, tolL, tolU, tolC] = isoutlier(acceptance, 'ThresholdFactor', 2); % Identify outliers
TF = acceptance < min(max(tolL, 0.1), tolC); % Thresholding low acceptance rates
badchains = find(TF); % Identify bad chains
% Plotting the acceptance rate thresholds and bad chains
yline(tolL, 'b--');
yline(tolU, 'b--');
yline(tolC, 'r--');
uq_postProcessInversion(myBayesianAnalysis_Step, 'badChains', badchains);
scatter(badchains, acceptance(badchains), 'red', 'filled');
legend('All chains', 'Lower Bound', 'Upper Bound', 'Median');

% Reporting and Displaying Bayesian Analysis Results
uq_print(myBayesianAnalysis_Step); % Print out the analysis results
uq_display(myBayesianAnalysis_Step); % Display results (prior and posterior distributions)

% Extracting optimal parameters from the Bayesian analysis
optpar_tmp = myBayesianAnalysis_Step.Results.PostProc.PointEstimate.X{1, 1};
optpar_tmp = 10.^optpar_tmp(1:4); % Transforming the parameters back from log scale

% Saving the Bayesian analysis results
filename_mat = sprintf('SaveResults_IUQ_GasAdsorption_%s.mat',ModelName);
save(filename_mat); % Save the current results

% Forward Modeling using Inferred Parameters
for i = 1:NumCurve
    Dm = optpar_tmp(3);
    Sc = mu/(rho_g*Dm);
    Re = us(i)*dp*rho_g/mu;
    Dz = Dm * (20+ 0.5*Sc*Re)/epsilon;
    
    kg = Dm/dp*(2.0+ 1.8 * Re^0.5 * Sc^(1/3));
    Rp = dp/2;
    De = optpar_tmp(4);
    q0star =Isotherm_Langmuir(cFeed,optpar_tmp(1),optpar_tmp(2));
    inv_K = (Rp*rho_b(i)*q0star)./(3*kg*cFeed*epsilon) + (Rp^2*rho_b(i)*q0star)./(15*De*cFeed*epsilon);
    K_G = 1/inv_K;
    
    par = [K_G,optpar_tmp(1),optpar_tmp(2),Dz]; %disp(par)
    
    optpar_disp(i,:) = par;
    
    Ngrid = ngrid(i);
    
    [t,c,~]= Model_GasAdsorption(10,Ngrid,dz,log10(par),epsilon,cFeed,200*60,ModelNum,us(i),rho(i));
    
    c_breakthrough=(c(Ngrid,:)+c(Ngrid+1,:))./2;
    c_est = c_breakthrough/cFeed;
    
    figure;
    filename = sprintf('ASZM-TEDA+GB(%.1fsccm,%dcm).csv',Q_read(i),L_read(i));
    exp_inf= csvread(filename,1,0);
    t_exp_tmp = exp_inf(:,1);
    c_exp_tmp = exp_inf(:,2);
    hold on
    plot(t_exp_tmp,c_exp_tmp,'o')
    hold on
    plot(t/60,c_est)
    xlabel('t [min]')
    ylabel('c/c0')
    ylim([0 1])
end

% Saving the results of forward modeling
filename_csv = sprintf('IUQ_Results_OptimalParameter_%s.csv',ModelName);
dlmwrite(filename_csv, optpar_disp, 'delimiter', ',', '-append'); % Append results to a CSV file