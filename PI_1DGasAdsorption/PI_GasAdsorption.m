%% General Code Description of  PI_GasAdsorption.m and Author Introduction
%======================================================================================================================
%> @ Code Description:
%> @ File       PI_GasAdsorption.m
%> @ Brief      An advanced parameter identification module based on Bayesian inference for gas adsorption processes.
%> @ Details    PI_GasAdsorption is the main program for parameter identification in gas adsorption model.
%>              This module includes loading of experimental data and setting up the Bayesian inference framework.
%>
%>              The parameter identification process in gas adsorption models via Bayesian inference unfolds as follows:
%>              (1) Parameter Establishment: Starts by defining a range of parameters, grounded in the physical attributes of gases and adsorbents.
%>              (2) Bayesian Inference Optimization: Adjusts the Bayesian inference settings in alignment with the unique properties of the materials and experimental conditions.
%>                  This includes optimizing the number of steps, samples, and the calibration of model predictions with actual measurements for enhanced accuracy.
%>              (3) Iterative Refinement and Analysis: This meticulous approach ensures that each iteration of Bayesian inference is both efficient and effective.
%>                  The process culminates in a comprehensive analysis, where the estimated parameters from the model are compared with the actual experimental data, often visualized through detailed graphs.

%>
%> @ Author Introduction
%>              - Yesol Hyun (yesol2@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University
%>              - Geunwoo Oh (gwoh@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University
%>              - Jung-Il Choi (jic@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University

%> @ Date        November 2023
%> @ Version     1.0
%> @ Cpyright   Copyright (c) 2023-2023 Yesol Hyun, Geunwoo Oh, and Jung-Il Choi, Yonsei University, All rights reserved.
%> @ License     This project is released under the terms of the MIT License (see LICENSE file).
%======================================================================================================================

%> @brief       Module for Bayesian parameter identification in gas adsorption models.
%> @details     This module facilitates parameter identification for gas adsorption processes using Bayesian inference.
%>              The process involves comparing experimental data with simulation results from the gas adsorption model.
%>              The steps include setting initial values and options for Bayesian inference, loading experimental data,
%>              numerically solving the model with estimated parameters, and optimizing these parameters through Bayesian analysis.

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
    idx2 = find(c_data > 0 & c_data <= 0.95, 1, 'last');
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
ModelOpts1.mHandle = @(par) uq_PhysicalAdsorption(par, t_exp, c_exp, dt, ngrid, dz, cFeed, epsilon, tf, v,ModelNum,rho, rho_b, rho_g, dp, mu, NumCurve, max_dcdt, max_tbr, Ndata);
ModelOpts1.isVectorized = true;
myForwardModel1 = uq_createModel(ModelOpts1);
BayesOpts.ForwardModel = myForwardModel1;

% Defining prior distributions for Bayesian parameters
PriorOpts1.Marginals(1).Name = 'b1';
PriorOpts1.Marginals(1).Type = 'Uniform';
PriorOpts1.Marginals(1).Parameters = [1 5];

PriorOpts1.Marginals(2).Name = 'qm';
PriorOpts1.Marginals(2).Type = 'Uniform';
PriorOpts1.Marginals(2).Parameters = [-2 0];

PriorOpts1.Marginals(3).Name = 'Dm';
PriorOpts1.Marginals(3).Type = 'Uniform';
PriorOpts1.Marginals(3).Parameters = [-7 -5];

PriorOpts1.Marginals(4).Name = 'De';
PriorOpts1.Marginals(4).Type = 'Uniform';
PriorOpts1.Marginals(4).Parameters = [-7 -5];

% Creating the prior distribution as an input object
myPriorDist1 = uq_createInput(PriorOpts1);
BayesOpts.Prior = myPriorDist1;

% Setting up the discrepancy model for the Bayesian analysis
SigmaOpts1.Marginals.Name = 'sigma2';
SigmaOpts1.Marginals.Type = 'Uniform';
SigmaOpts1.Marginals.Parameters = [0 (0.05)^2];
mySigmaDist = uq_createInput(SigmaOpts1);
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
myBayesianAnalysis_Step1 = uq_createAnalysis(BayesOpts); % Create Bayesian analysis
uq_time = toc; % Record the time taken for the analysis

% Filtering MCMC chains with low acceptance rates
uq_display(myBayesianAnalysis_Step1, 'acceptance', true); % Display acceptance rates
acceptance = myBayesianAnalysis_Step1.Results.Acceptance; % Extract acceptance rates
[~, tolL, tolU, tolC] = isoutlier(acceptance, 'ThresholdFactor', 2); % Identify outliers
TF = acceptance < min(max(tolL, 0.1), tolC); % Thresholding low acceptance rates
badchains = find(TF); % Identify bad chains
% Plotting the acceptance rate thresholds and bad chains
yline(tolL, 'b--');
yline(tolU, 'b--');
yline(tolC, 'r--');
uq_postProcessInversion(myBayesianAnalysis_Step1, 'badChains', badchains);
scatter(badchains, acceptance(badchains), 'red', 'filled');
legend('All chains', 'Lower Bound', 'Upper Bound', 'Median');

% Reporting and Displaying Bayesian Analysis Results
uq_print(myBayesianAnalysis_Step1); % Print out the analysis results
uq_display(myBayesianAnalysis_Step1); % Display results (prior and posterior distributions)

% Extracting optimal parameters from the Bayesian analysis
optpar_tmp = myBayesianAnalysis_Step1.Results.PostProc.PointEstimate.X{1, 1};
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
    
    [t,c,~]= Adsorption_Model(10,Ngrid,dz,log10(par),epsilon,cFeed,200*60,ModelNum,us(i),rho(i));
    
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
%% IUQ Procedure Module for Parameter Estimation in Physical Adsorption
%======================================================================================================================
%> @details     This function, uq_PhysicalAdsorption, is designed to execute the Inverse Uncertainty Quantification (IUQ) procedure for parameter estimation in physical adsorption models.
%>              It processes multiple data sets to estimate adsorption model parameters using Bayesian analysis.
%>              The function includes:

%>              (1) Parameter Transformation and Looping Over Data Sets:
%>                  - Transforms parameters from logarithmic to linear scale.
%>                  - Iterates over multiple experimental data curves.
%>
%>              (2) Parameter Calculation for Each Curve:
%>                  - Calculates adsorption-related parameters like mass transfer coefficients.
%>                  - Determines the breakthrough curve using a forward adsorption model.
%>
%>              (3) Data Interpolation and Differentiation:
%>                  - Interpolates the estimated concentration data to reference time points.
%>                  - Calculates the rate of change of concentration.
%>
%>              (4) Error Handling and Output Data Compilation:
%>                  - Implements checks for unrealistic model outputs.
%>                  - Compiles the differential rate and time-to-breakthrough data for each parameter set.
%======================================================================================================================
function t_IUQ = uq_PhysicalAdsorption(X, tref_multi, cref_multi, dt, ngrid, dz, cFeed, epsilon, tf, us,num,rho, rho_b, rho_g, dp, mu, NumCurve, max_dcdt, max_tbr, Ndata)

Npar = size(X, 1); % Number of parameter sets
for j = 1:Npar
    % Transform parameters from log scale to actual scale
    par_tmp = 10.^X(j, 1:4);
    
    for i = 1:NumCurve
        c_ref = cref_multi(1+Ndata*(i-1):Ndata+Ndata*(i-1));
        t_ref = tref_multi(1+Ndata*(i-1):Ndata+Ndata*(i-1));
        
        Dm = par_tmp(3);
        Sc = mu/(rho_g*Dm);
        Re = us(i)*dp*rho_g/mu;
        Dz = Dm * (20+ 0.5*Sc*Re)/epsilon;
        
        kg = Dm/dp*(2.0+ 1.8 * Re^0.5 * Sc^(1/3));
        Rp = dp/2;
        De = par_tmp(4);
        q0star =Isotherm_Langmuir(cFeed,par_tmp(1),par_tmp(2));
        inv_K = (Rp*rho_b(i)*q0star)./(3*kg*cFeed*epsilon) + (Rp^2*rho_b(i)*q0star)./(15*De*cFeed*epsilon);
        K_G = 1/inv_K;
        par = [K_G,par_tmp(1),par_tmp(2),Dz]; %disp(par)
        
        Ngrid = ngrid(i);
        
        [t_est,c_tmp,~]= Adsorption_Model(dt,Ngrid,dz,log10(par),epsilon,cFeed,tf(i),num,us(i),rho(i));
        
        c_breakthrough=(c_tmp(Ngrid,:)+c_tmp(Ngrid+1,:))./2;
        c_est = c_breakthrough;
        c_interp = interp1(t_est,c_est,t_ref,'linear', 'extrap');
        
        dt = t_ref(2) - t_ref(1);
        dc_dt = diff(c_interp) / dt;
        
        if length(t_est)< 3 ||max(c_est)>cFeed || min(c_est)<0
            tmp_dcdt = 100*ones(1,length(c_ref)-1);
            tmp_tbr(i)= 100;
        else
            tmp_dcdt= dc_dt/max_dcdt(i);
            idx1=find(c_est> 0,1,'first');
            idx2=find(c_est> 0 & c_est< 0.95* cFeed,1,'last');
            interp_func = @(c) interp1(c_est(idx1:idx2)/cFeed, t_est(idx1:idx2)/60, c, 'linear', 'extrap');
            specific_c = 0.01;
            specific_t = interp_func(specific_c);
            tmp_tbr(i)= specific_t;
        end
        
        if i==1
            tmp_stack =tmp_dcdt;
        else
            tmp_stack = [tmp_stack,tmp_dcdt];
        end
    end
    
    % Compile differential rate and time-to-breakthrough data
    t_IUQ(j, :) = [((Ndata - 1) / Ndata) * tmp_stack, (1 / Ndata) * tmp_tbr / max_tbr];
end
end
%% Forward Modeling for Gas Physical Adsorption
%======================================================================================================================
%> @details     This function, Adsorption_Model, performs forward modeling for gas physical adsorption. It simulates the adsorption process using various isotherm models and solves the system using an implicit solver with adaptive time stepping. Key steps include:

%>              (1) Model Initialization:
%>                  - Transforms model parameters from logarithmic to linear scale.
%>                  - Sets up initial conditions and parameters based on the chosen isotherm model.
%>
%>              (2) Time-stepping and Implicit Solving:
%>                  - Iterates through time steps, updating concentrations and adsorbed quantities.
%>                  - Applies a finite volume method with implicit solving for the nonlinear model.
%>
%>              (3) Boundary Conditions and Matrix Assembly:
%>                  - Implements boundary conditions and assembles the system matrix for solving.
%>
%>              (4) Newton Iteration and Convergence:
%>                  - Performs Newton iteration for nonlinear system convergence.
%>                  - Checks and records convergence metrics.
%>
%>              (5) Results Compilation and Breakthrough Detection:
%>                  - Compiles concentration and adsorption data at each time step.
%>                  - Detects breakthrough times based on concentration thresholds.
%======================================================================================================================
% Function definition for the adsorption model
function [t, c, q] = Adsorption_Model(dt, ngrid, dz, par, epsilon, cFeed, tf,us, rho)
% Define model parameter 
par= 10.^par;

if num==1 %num 1= Langmuir
    k_LDF= par(1);
    b=par(2);
    qm=par(3);
    D= par(4);
elseif num==2 %num 2= Linear
    k_LDF=par(1);
    K_H=par(2);
    D= par(3);
elseif num==3 %num 3= Toth
    k_LDF=par(1);
    b=par(2);
    qm=par(3);
    n_Toth= par(4);
    D= par(5);
end

% Time stepping and implicit solving
t_tmp = 0:dt:tf;
ntime = length(t_tmp);
idx= zeros(2,1);
t= zeros(ntime,1);
c = zeros(ngrid+1,ntime); 
c(1,1) = cFeed;
q = zeros(ngrid+1,ntime); 
tol_save = zeros(ntime,100);

alpha= rho*(1-epsilon)/epsilon;
ui= us/epsilon; % interstitial velocity
beta = ui*dt/(2*dz);
gamma = D*dt/(2*dz*dz); 
lambda = k_LDF*dt/2;
Aim = -(beta + gamma); 
Aip = - gamma; 


for nt = 2:ntime
    t(nt)=t(nt-1)+dt;
    t_cn= (t(nt)+t(nt-1))/2;
    
    % Preallocate  and initialize necessary vector
    n = ngrid-1;
    Amatrix= zeros(n);
    delta_old= zeros(n,1);
    F= zeros(n,1);
    
    c_old = c(:,nt-1);
    q_old= q(:,nt-1);
    
    if nt==2
        if num==1
            g_old=Isotherm_Langmuir(c_old,b,qm);
        elseif num==2
            g_old=Isotherm_Linear(c_old,K_H);
        elseif num==3
            g_old=Isotherm_Toth(c_old,b,qm,n_Toth); 
        end
    end
    
    % Newton iteration for solving the nonlinear governing equations 
    for iter = 1:5
        if iter==1
            c_tmp =  c_old;
        end
        
        for nz=2:ngrid
            if num==1
                g_tmp=Isotherm_Langmuir(c_tmp(nz),b,qm);
                gprime_tmp=Deriv_Langmuir(c_tmp(nz),b,qm);
            elseif num==2
                g_tmp=Isotherm_Linear(c_tmp(nz),K_H);
                gprime_tmp=Deriv_Linear(K_H);
            elseif num==3
                g_tmp=Isotherm_Toth(c_tmp(nz),b,qm,n_Toth);
                gprime_tmp=Deriv_Toth(c_tmp(nz),b,qm,n_Toth);
            end
            
            F_conv= beta*(c_tmp(nz)-c_tmp(nz-1)+c_old(nz)-c_old(nz-1));
            F_diff= gamma*(c_tmp(nz+1)-2*c_tmp(nz)+c_tmp(nz-1)+c_old(nz+1)-2*c_old(nz)+c_old(nz-1));
            q_tmp= (lambda*(g_tmp+ g_old(nz))+(1-lambda)*q_old(nz))/(1+lambda);
            F_kinetic= alpha*lambda*(g_tmp+g_old(nz)- q_tmp-q_old(nz));
            
            F(nz-1)= c_tmp(nz)-c_old(nz) + F_conv - F_diff + F_kinetic;
            
            Aic = 1 + beta + 2*gamma + alpha*lambda*gprime_tmp/(1+lambda);
            Amatrix(nz-1,nz-1) = Aic;
            
            if(nz-1 ~= 1)
                Amatrix(nz-1,nz-2) = Aim;
                Amatrix(nz-2,nz-1) = Aip;
            end
            
            % Boundary treatment 
            Amatrix(1,1) = Amatrix(1,1) + gamma;
        end      
        
        RHS = -F;
        
        if iter==1
            max_Fold = max(F);
        end
        
        % Solving the system using Thomas algorithm
        delta = thomas(Amatrix, RHS);
        
        % Updating concentration and checking convergence
        c_tmp(2:ngrid)= c_tmp(2:ngrid)+delta;
        c_tmp(1) = 2*cFeed  - c_tmp(2);
        c_tmp(ngrid+1) = c_tmp(ngrid);
        tol= sqrt(sum((delta-delta_old).^2)/length(delta));
        tol_save(nt,iter)= tol;
        delta_old= delta;
        
        max_F= max(F);
        relF= abs(max_F/max_Fold);
        
        % Break from iteration upon convergence
        if (relF <= 1e-10 || iter >= 5)
            c_fin = c_tmp;
            break;
        end
    end
    
    % Update concentration and adsorption quantities for current time step
    c(:, nt) = c_fin;

    if num==1
        g_new=Isotherm_Langmuir(c(:,nt),b,qm);
    elseif num==2
        g_new=Isotherm_Linear(c(:,nt),K_H);
    elseif num==3
        g_new=Isotherm_Toth(c(:,nt),b,qm,n_Toth);
    end
    
    q(:,nt) = (lambda*(g_new+g_old) + (1-lambda)*q_old)/(1+lambda); % calculating new q
    g_old = g_new;
       
    % Detection of breakthrough times
    if (c(ngrid, nt) / cFeed >= 0.001 && c_old(ngrid) / cFeed <= 0.001)
        idx(1) = nt + 1;
    end
    if (t_cn > tf)
        idx(2) = nt;
    end
end
end

%% Supporting Functions for Adsorption Modeling
%======================================================================================================================
%> @details     This section includes supporting functions used in the adsorption modeling process, specifically for solving linear systems and calculating isotherm properties. The functions include:

%>              (1) Thomas Algorithm Function:
%>                  - Efficiently solves tridiagonal linear systems Ax=d using the Thomas algorithm.
%>                  
%>              (2) Langmuir Isotherm Function:
%>                  - Calculates the adsorption quantity based on the Langmuir isotherm model.
%>                  
%>              (3) Derivative of Langmuir Isotherm Function:
%>                  - Computes the derivative of the Langmuir isotherm model with respect to concentration.
%======================================================================================================================
%% Thomas Algorithm Function
function x = thomas(A, d)
% Solves a tridiagonal linear system Ax=d
n = length(d); % Determine the size of the system
% Preallocate vectors for the sub-diagonal, diagonal, and super-diagonal
a = zeros(n-1, 1);
b = zeros(n, 1);
c = zeros(n-1, 1);
x = zeros(n, 1);

% Initialize the diagonal
b(1) = A(1, 1);

% Forward elimination process
for i = 2:n
    % Extract elements from A and perform elimination
    a(i-1) = A(i, i-1);
    b(i) = A(i, i);
    c(i-1) = A(i-1, i);
    w = a(i-1) / b(i-1);
    b(i) = b(i) - w * c(i-1);
    d(i) = d(i) - w * d(i-1);
end

% Backward substitution process
x(n) = d(n) / b(n);
for i = (n-1):-1:1
    x(i) = (d(i) - c(i) * x(i+1)) / b(i);
end
end
%% Langmuir Isotherm Function
function g = Isotherm_Langmuir(c, b, qm)
% Calculates the Langmuir isotherm
g = qm * (b * c) ./ (1 + b * c);
end
%% Derivative of Langmuir Isotherm Function
function gprime = Deriv_Langmuir(c, b, qm)
% Calculates the derivative of the Langmuir isotherm
gprime = (b * qm) ./ (b * c + 1) - (b^2 * c * qm) ./ (b * c + 1).^2;
end
