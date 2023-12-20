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
function [t, c, q] = Model_GasAdsorption(dt, ngrid, dz, par, epsilon, cFeed, tf,num,us, rho)
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
%% Langmuir isotherm
function g =Isotherm_Langmuir(c,b,qm)
g= qm*(b*c)./(1+b*c);
end
%% Freudlich
function g =Isotherm_Freudlich(c,b,qm)
c= max(c,0);
g=qm*real(c.^(1/b));
end
%% Linear
function g =Isotherm_Linear(c,qm)
g=qm*c;
end
%% Toth
function g =Isotherm_Toth(c,b,qm,t)
c= max(c,0);
g= qm*(b*c)./((1+(b*c).^t).^(1/t));
end
%% Derivative of Langmuir isotherm
function gprime =Deriv_Langmuir(c,b,qm)
gprime= (b*qm)./(b*c + 1) - (b^2*c*qm)./(b*c + 1).^2 ;
end
%% Derivative of Freudlich isotherm
function gprime =Deriv_Freudlich(c,b,qm)
c= max(c,0);
gprime=(c.^(1/b - 1)*qm)/b;
end
%% Derivative of Linear isotherm
function gprime =Deriv_Linear(qm)
gprime=qm;
end
%% Derivative of Toth isotherm
function gprime =Deriv_Toth(c,b,qm,t)
c= max(c,0);
gprime= (b*qm)/((b*c).^t + 1).^(1/t) - (b^2.*c*qm*(b*c).^(t - 1))/((b.*c).^t + 1).^(1/t + 1);
end