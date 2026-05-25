function [Tbase_list] = latin_hypercube_sampling(sampling_number)
    %============Inputs============
    %**System Paremeters**
    TDP = 125.0;           % Thermal Design Power, W
    NumFins = 5;           % Number of heat sink fins
    L = 8/100;             % Total length of domain (CPU + heat sink), m
    L_CPU = 3/100;         % Length of CPU region, m
    thick = 0.3/100;       % Thickness of heat sink fin, m

    %**Material Parameters**
    k_Metal = 160.;               % Thermal conductivity of Metal, W/(m*K)
    density_Metal = 2700.;        % Density of metal, kg/m^3
    specific_heat_Metal = 895.0;  % Specific heat of metal, J/(kg*K)
    alpha_Metal = k_Metal/(density_Metal*specific_heat_Metal);   % Thermal diffusivity of Metal, m^2/s
    %====NOTE: k_CPU is the uncertain input variable
    %===========================================================
    k_CPU =0.4;                   % Thermal conductivity of CPU, W/(m*K)
    %===========================================================
    density_CPU = 2000.;          % Density of CPU, kg/m^3
    specific_heat_CPU = 750.;     % Specific heat of CPU, J/(kg*K)


    %**Model Parameters (averaged in flow direction)**
    hbar = 15.0;           % Averaged convective cooling coeff., W/(m^2*K)
    Tbar = 330.0;          % Averaged air temperature, K

    %**Discretization Parameters**
    imax = 33;             % Number of nodes in x (normal to flow direction)
    itermax = 1e7;         % Maximum allowable number of iterations
    convtol = 1.e-8;       % Iterative convergence tolerance (relative to fifth iteration)
    
    i=1;
    sigma=[0.0049, 0; 0, 0.0025];
    mu=[0.4, 0.3];
    values=[];
    lhs_=lhsnorm(mu, sigma, sampling_number);
    while i<= sampling_number
        temp = lhs_(i, :);
        random_kcpu=temp(1)
        random_thick=temp(2)
        random_thick=random_thick/100;
        T(1:imax) = 350.0;     % Temperature of heat sink (initialized to 350 K), K
        alpha_CPU = random_kcpu/(density_CPU*specific_heat_CPU); % Thermal diffusivity of CPU, m^2/s
        [Tbase, T, x, L2conv, history] = heatcondsolve(T,TDP,NumFins,L,random_thick,L_CPU,k_Metal,alpha_Metal,random_kcpu,alpha_CPU,hbar,Tbar,imax,itermax,convtol);
        values=[values, Tbase-273.15];
        i = i + 1;
    end
    Tbase_list=sort(values);
 end

