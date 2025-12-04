clc
close all
clear variables
tic                    % Timer Start

%% To Do:


%% Test Case 1 (YouTube @giuseppe):
% R1_vec = [0;0;0];             % Placeholder
% R2_vec = [0;0;0];             % Placeholder
% r1 = 1.006;                   % AU
% r2 = 0.727;                   % AU
% theta = 264.8 + 9.1;          % deg
% 
% ToF_input = 158;              % day
% Long = true;                  % Long Path T/F
% N = 0;                        % Number of Revolutions
% mu = 2.95912208339 * 10 ^ -4; % AU^3/day^2 for the Sun
% 
% max_iterations = 100;

% Expected: Semi-major axis = 0.7655 AU, alpha = 3.5544, beta = -0.8604
% Actual: Semi-major axis = 0.7655 AU, alpha = 3.5544, beta = -0.8604

%% Test Case 2 (Prussing Long Fig. 3):
% R1_vec = [0;0;0];             % Placeholder
% R2_vec = [0;0;0];             % Placeholder
% r1 = 1.0;                     % Length Unit
% r2 = 1.524;                   % Length Unit
% theta = 107;                  % deg
% 
% ToF_input = 1;                % Canonical Time Unit?
% Long = true;                  % Long Path T/F
% N = 0;                        % Number of Revolutions
% mu = 4*pi*pi;                 % Canonical Unit?
% 
% max_iterations = 100;

% Expected: Semi-major axis == 1.25 Length Unit
% Actual: Semi-major axis = 1.2752 Length Unit

%% Test Case 3 (Prussing Short Fig. 3):
% R1_vec = [0;0;0];             % Placeholder
% R2_vec = [0;0;0];             % Placeholder
% r1 = 1.0;                     % Length Unit
% r2 = 1.524;                   % Length Unit
% theta = 107;                  % deg
% 
% ToF_input = 0.5;              % Canonical Time Unit?
% Long = false;                 % Long Path T/F
% N = 0;                        % Number of Revolutions
% mu = 4*pi*pi;                 % Canonical Unit?
% 
% max_iterations = 100;

% Expected: Semi-major axis == 1.17 Length Unit
% Actual: Semi-major axis = 1.1731 Length Unit

%% Test Case 4 (Shen Short Fig. 2):
% R1_vec = [0;0;0];             % Placeholder
% R2_vec = [0;0;0];             % Placeholder
% r1 = 1.0;                     % Length Unit
% r2 = 2;                       % Length Unit
% theta = 60;                   % deg
% 
% ToF_input = 0.5;              % Canonical Time Unit?
% Long = false;                 % Long Path T/F
% N = 0;                        % Number of Revolutions
% mu = 4*pi*pi;                 % Canonical Unit?
% 
% max_iterations = 100;

% Expected: Semi-major axis == 1.21 Length Unit
% Actual: Semi-major axis = 1.2088 Length Unit

%% Test Case 5 (Shen Long Fig. 2):
% R1_vec = [0;0;0];             % Placeholder
% R2_vec = [0;0;0];             % Placeholder
% r1 = 1.0;                     % Length Unit
% r2 = 2;                       % Length Unit
% theta = 60;                   % deg
% 
% ToF_input = 0.75;             % Canonical Time Unit?
% Long = true;                  % Long Path T/F
% N = 0;                        % Number of Revolutions
% mu = 4*pi*pi;                 % Canonical Unit?
% 
% max_iterations = 100;

% Expected: Semi-major axis == 1.21 Length Unit
% Actual: Semi-major axis = 1.2099 Length Unit

%% Test Case 6 (Shen Long Multi Rev Fig. 2):
% R1_vec = [0;0;0];             % Placeholder
% R2_vec = [0;0;0];             % Placeholder
% r1 = 1.0;                     % Length Unit
% r2 = 2;                       % Length Unit
% theta = 60;                   % deg
% 
% ToF_input = 6.25;             % Canonical Time Unit?
% Long = true;                  % Long Path T/F
% N = 2;                        % Number of Revolutions
% mu = 4*pi*pi;                 % Canonical Unit?
% 
% max_iterations = 100;

% Expected: Semi-major axis == 1.7 Length Unit
% Actual: Semi-major axis = 1.702 Length Unit

%% Test Case 7 (Shen Short Multi Rev Fig. 2):
% R1_vec = [0;0;0];             % Placeholder
% R2_vec = [0;0;0];             % Placeholder
% r1 = 1.0;                     % Length Unit
% r2 = 2;                       % Length Unit
% theta = 60;                   % deg
% 
% ToF_input = 7.6;              % Canonical Time Unit?
% Long = false;                 % Long Path T/F
% N = 3;                        % Number of Revolutions
% mu = 4*pi*pi;                 % Canonical Unit?
% 
% max_iterations = 100;

% Expected: Semi-major axis == 1.8 Length Unit
% Actual: Semi-major axis = 1.8056 Length Unit

%% Test Case 8 (Superior Prograde Ex. 1):
% R1_vec = [22592.145603; -1599.915239; -19783.950506];   % km
% R2_vec = [1922.067697; 4054.157051; -8925.727465];      % km
% 
% r1 = norm(R1_vec);                                      % km
% r2 = norm(R2_vec);                                      % km
% theta = acosd(dot(R1_vec,R2_vec)/(r1*r2));              % deg
% 
% ToF_input = 36000;                                      % sec
% Long = true;                                            % Long Path T/F
% N = 0;                                                  % Number of Revolutions
% mu = 3.986004418*10^5;                                  % Earth Grav Parameter km3 s-2
% 
% max_iterations = 50;

% Expected: V1 == [2.000652697; 0.387688615; -2.666947760],
% V2 == [-3.79246619; -1.77707641; 6.856814395]
% Actual: V1 = [2.00065269702563; 0.387688615292766; -2.66694775975572],
% V2 = [-3.79246618851045; -1.77707640626944; 6.85681439477726]

%% Test Case 9 (Superior Retrograde Ex. 1):
% R1_vec = [22592.145603; -1599.915239; -19783.950506];   % km
% R2_vec = [1922.067697; 4054.157051; -8925.727465];      % km
% 
% r1 = norm(R1_vec);                                      % km
% r2 = norm(R2_vec);                                      % km
% theta = 360-acosd(dot(R1_vec,R2_vec)/(r1*r2));          % deg
% 
% ToF_input = 36000;                                      % sec
% Long = true;                                            % Long Path T/F
% N = 0;                                                  % Number of Revolutions
% mu = 3.986004418*10^5;                                  % Earth Grav Parameter km3 s-2
% 
% max_iterations = 50;

% Expected: V1 == [2.96616042; -1.27577231; -0.75545632],
% V2 == [5.84375455; -0.20047673; -5.48615883]
% Actual: V1 = [2.96616041748762; -1.27577231102832; -0.755456316568433],
% V2 = [5.84375454679398; -0.200476733558268; -5.48615882868479]

%% Test Case 10 (Superior Ex. 2):
% R1_vec = [7231.58074563487; 218.02523761425; 11.79251215952];   % km
% R2_vec = [7357.06485698842; 253.55724281562; 38.81222241557];   % km
% 
% r1 = norm(R1_vec);                                              % km
% r2 = norm(R2_vec);                                              % km
% theta = acosd(dot(R1_vec,R2_vec)/(r1*r2));                      % deg
% 
% ToF_input = 12300;                                              % sec
% Long = true;                                                    % Long Path T/F
% N = 0;                                                          % Number of Revolutions
% mu = 3.986004418*10^5;                                          % Earth Grav Parameter km3 s-2
% 
% max_iterations = 100;                                           % Important to keep high, insane sma

% As Expected, change sma_max to 1e+07

%% Test Case 11 (Real Mars 26-Jun-2035 to 04-Jan-2036):
% R1_vec = [1.12826e+07; -1.51646e+08; 11596];        % km (Earth Dept.)
% R2_vec = [1.11252e+08 1.94131e+08 1.34136e+06];   % km (Mars Arrival)
% r1 = norm(R1_vec);                                  % km
% r2 = norm(R2_vec);                                  % km
% theta = acosd(dot(R1_vec,R2_vec)/(r1*r2));          % deg
% 
% ToF_input = 192*24*60*60;                           % sec
% Long = false;                                       % Long Path T/F
% N = 0;                                              % Number of Revolutions
% mu = 1.327124e11;                                   % km3 s-2
% 
% max_iterations = 100;
% 
% sma_Earth = 149.6e6; % km
% sma_Mars = 228e6; % km

%% Test 12 (Mars N=0 Short Optimal):
R1_vec = [9.37851e+07 1.14552e+08 -10416.8];        % km (Earth Dept.)
R2_vec = [1.55349e+08 1.53576e+08 -589769];   % km (Mars Arrival)
r1 = norm(R1_vec);                                  % km
r2 = norm(R2_vec);                                  % km
theta = acosd(dot(R1_vec,R2_vec)/(r1*r2));          % deg

ToF_input = 1080*24*60*60;                           % sec
Long = true;                                       % Long Path T/F
N = 1;                                              % Number of Revolutions
mu = 1.327124e11;                                   % km3 s-2

max_iterations = 100;

sma_Earth = 1.49598e+08; % km
sma_Mars = 2.27933e+08; % km




%% Calcs
if theta > 180
    beta_sign = -1;
else
    beta_sign = 1;
end

chord = sqrt(r1*r1 + r2*r2 - 2*r1*r2*cosd(theta)); % Length
semiperimeter = (r1+r2+chord) * 0.5;               % Length

sqrt_mu = sqrt(mu);

semi_major_axis_max = semiperimeter*10000;               % Length (should be infinity)
semi_major_axis_min = semiperimeter/2;             % Length



%% Loop

disp('Beginning binary search')
for i = 1:max_iterations

    % Mid Calcs
    semi_major_axis_mid = (semi_major_axis_min + semi_major_axis_max) * 0.5;                  % Length
    
    if Long == true
        alpha_mid = 2*pi-2*asin(sqrt(semiperimeter/(semi_major_axis_mid*2)));                 % rad
    else
        alpha_mid = 2*asin(sqrt(semiperimeter/(semi_major_axis_mid*2)));                      % rad
    end
    
    beta_mid = beta_sign*2*asin(sqrt((semiperimeter-chord)/(semi_major_axis_mid*2)));         % rad

    ToF_mid = (sqrt(semi_major_axis_mid^3)*(2*N*pi+alpha_mid-beta_mid-(sin(alpha_mid)-sin(beta_mid))))/sqrt_mu; % Time


    % Max Calcs
    if Long == true
        alpha_max = 2*pi-2*asin(sqrt(semiperimeter/(semi_major_axis_max*2)));                 % rad
    else
        alpha_max = 2*asin(sqrt(semiperimeter/(semi_major_axis_max*2)));                      % rad
    end
    
    beta_max = beta_sign*2*asin(sqrt((semiperimeter-chord)/(semi_major_axis_max*2)));         % rad

    ToF_max = (sqrt(semi_major_axis_max^3)*(2*N*pi+alpha_max-beta_max-(sin(alpha_max)-sin(beta_max))))/sqrt_mu; % Time


    % Min Calcs
    if Long == true
        alpha_min = 2*pi-2*asin(sqrt(semiperimeter/(semi_major_axis_min*2)));                 % rad
    else
        alpha_min = 2*asin(sqrt(semiperimeter/(semi_major_axis_min*2)));                      % rad
    end
    
    beta_min = beta_sign*2*asin(sqrt((semiperimeter-chord)/(semi_major_axis_min*2)));         % rad

    ToF_min = (sqrt(semi_major_axis_min^3)*(2*N*pi+alpha_min-beta_min-(sin(alpha_min)-sin(beta_min))))/sqrt_mu; % Time


    % Comparison
    ToF_List_small = [ToF_min, ToF_mid];
    ToF_List_big = [ToF_mid, ToF_max];

    if isbetween(ToF_input,min(ToF_List_small),max(ToF_List_small))
        semi_major_axis_max = semi_major_axis_mid;

    elseif isbetween(ToF_input,min(ToF_List_big),max(ToF_List_big))
        semi_major_axis_min = semi_major_axis_mid;
    else
        error('ToF Range Error, check sma min and max')
    end
end

%% Eccentricity Calc

eccentricity = sqrt(1 - (((4 * (semiperimeter - r1) * (semiperimeter - r2)) / (chord*chord)) * sin((alpha_mid + beta_mid) / 2) * sin((alpha_mid + beta_mid) / 2)));

%% Velocity Calc

if norm(R1_vec) ~= 0
    A = sqrt(mu/(4*semi_major_axis_mid))*cot(alpha_mid/2);
    B = sqrt(mu/(4*semi_major_axis_mid))*cot(beta_mid/2);
    
    U1_vec = R1_vec./r1; 
    U2_vec = R2_vec./r2; 
    Uc_vec = (R2_vec - R1_vec)./chord;
    
    V1_vec = (B+A)*Uc_vec + (B-A)*U1_vec;
    V1_norm = norm(V1_vec);
    V2_vec = (B+A)*Uc_vec - (B-A)*U2_vec;
    V2_norm = norm(V2_vec);
end

if sma_Earth > 0
    V1_pretransfer = sqrt(mu*((2/r1)-(1/sma_Earth)));
    V2_pretransfer = sqrt(mu*((2/r2)-(1/sma_Mars)));
end

%% Outputs

disp('All done!')
disp(['Converged sma = ', num2str(semi_major_axis_mid), ', ToF = ', num2str(ToF_mid), ', Number of Revolutions = ', num2str(N)])

if norm(R1_vec) ~= 0
    disp(['Velocity at Position 1 = ',mat2str(V1_vec), ', V1 mag = ',num2str(V1_norm), ', V1_Earth_mag = ',num2str(V1_pretransfer)])
    disp(['Velocity at Position 2 = ',mat2str(V2_vec), ', V2 mag = ',num2str(V2_norm), ', V2_Mars_mag = ',num2str(V2_pretransfer)])
end

disp(['Eccentricity = ', num2str(eccentricity), ', Alpha = ', num2str(alpha_mid), ' rad, Beta = ', num2str(beta_mid), ' rad, Theta = ', num2str(theta), ' deg'])
disp(['Runtime = ',num2str(toc), ' sec for ', num2str(i), ' iterations'])                     % Timer Stop