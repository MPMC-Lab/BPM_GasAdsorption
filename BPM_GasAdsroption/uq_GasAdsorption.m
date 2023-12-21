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
function t_IUQ = uq_GasAdsorption(X, tref_multi, cref_multi, dt, ngrid, dz, cFeed, epsilon, tf, us,num,rho, rho_b, rho_g, dp, mu, NumCurve, max_dcdt, max_tbr, Ndata,Max_C)

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
        
        [t_est,c_tmp,~]= Model_GasAdsorption(dt,Ngrid,dz,log10(par),epsilon,cFeed,tf(i),num,us(i),rho(i));
        
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
            idx2=find(c_est> 0 & c_est< Max_C* cFeed,1,'last');
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