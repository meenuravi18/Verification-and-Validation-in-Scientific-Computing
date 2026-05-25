function [Tbase, T, x, L2conv, history] = heatcondsolve(T,TDP,NumFins,L,thick,L_CPU,k_Metal,alpha_Metal,k_CPU,alpha_CPU,hbar,Tbar,imax,itermax,convtol)
% Solves the 1D heat equation with a specified energy flux to approximate 
% the base temperature of a heat sink with fins in a laminar fluid medium.
%
% [Tbase, T, x, L2conv, history] =
% heatcondsolve(T,TDP,NumFins,L,thick,L_CPU,k_Metal,alpha_Metal,k_CPU,alpha_CPU,hbar,Tbar,imax,itermax,convtol)
% 
%
% Inputs:
%                     T : Initial Temperature (matrix [1 x imax] or [imax x 1] )
%                   TDP : Thermal Design Power, W 
%               NumFins : Number of heat sink fins 
%                     L : Length of fins, m 
%                 thick : Thickness of heat sink fin, m 
%                 L_CPU : Lenth of CPU region, m 
%               k_Metal : Thermal conductivity of Metal, W/(m*K)
%           alpha_Metal : Thermal diffusivity of Metal, m^2/s
%                 k_CPU : Thermal conductivity of Metal, W/(m*K)
%             alpha_CPU : Thermal diffusivity of Metal, m^2/s
%                  hbar : Averaged convective cooling coeff., W/(m^2*K)
%                  Tbar : Averaged air temperature, K
%                  imax : number of nodes in x
%               itermax : maximum number of iterations
%               convtol : Normalized convergence tolerance
% 
% Outputs: 
%                 Tbase : Temperature at base of heat sink
%                     T : Final temperature
%                     x : x node locations
%                L2conv : Final normalized residual
%                istory : Iterative residual history  of every 100th iteration [iter #, residual]


%% Initialize variables
history = [0,0];        % Iterative residual history
L2normalize=1;          % Normalizing value for iterative residuals
insul = zeros(size(T)); % Insulation blanking array
k=insul;                % Thermal conductivity
alpha=insul;            % Thermal diffusivity

%% Error checking
if any(size(T)~= [1,imax]) && any(size(T)~=[imax,1])
    error('Temperature vector must match specified number of nodes!')
end

%% Set geometry
x=linspace(0,L,imax);
dx = L/(imax-1);

%% Set thermal conductivity and thermal diffusivity arrays
for i = 1:imax
    k(i) = k_CPU;
    alpha(i) = alpha_CPU;
    insul(i) = 0.;
    if x(i) > L_CPU 
        insul(i) = 1.0;
        k(i) = k_Metal;
        alpha(i) = alpha_Metal;
    end 
    if abs(x(i) - L_CPU)<1.e-6
        insul(i) = 0.5;
        k(i) = 0.5*(k_Metal + k_CPU);
        alpha(i) = 0.5*(alpha_Metal + alpha_CPU);
    end
end


%% Set time step array (local time stepping)
dt = 0.4*dx^2./alpha;


%% Iterative solution loop
i=2:imax-1;
src_constant = 2*alpha.*insul.*hbar./(k.*thick);
hist_cnt = 1;
for k = 1:itermax

    % Compute Right Hand Side
    d2Tdx2 = ( T(i+1) - 2*T(i) + T(i-1) )/(dx^2);
    src = src_constant(i).*(T(i)-Tbar);
    RHS(i) = alpha(i).*d2Tdx2 - src;
    L2conv = sqrt( sum(RHS(i).^2)/(imax-2) );
   
    % Update Temperatures in Interior
    T(i) = T(i) + dt(i).*RHS(i);

    % Update Boundary Conditions
    T(1) = ( 4*T(2) - T(3) + 2.*dx*TDP / ( k_CPU*NumFins ) ) / 3.0;
    
    psi = 2*dx*hbar*( dx/thick + 1 )/k_Metal;
    T(imax) = ( 4*T(imax-1) - T(imax-2) + psi*Tbar ) / ( 3.0 + psi );
    
    % Normalize residuals
    if (k==5); L2normalize=L2conv; end
    L2conv=L2conv/L2normalize;
    
    % Store iterative residuals
    if mod(k,100)==0; 
        history(hist_cnt,:) = [k,L2conv]; 
        hist_cnt=hist_cnt+1; 
    end

    % Check for iterative convergence
    if (L2conv < convtol) && (k>2); break ; end
        
end

%% Convergence checking
if (L2conv > convtol)
    fprintf('Warning! Solution not Converged. Imax = %4.0f, Residual = %5.3e\n',imax, L2conv)
end

%% Output
Tbase = T(1);

end
