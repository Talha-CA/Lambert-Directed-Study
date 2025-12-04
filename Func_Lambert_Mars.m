% Lambert Solver, rewritten as a function, requires R vectors and given N

% Inputs                : Units,                                Description
% R1_vec                : any Length units,                     Starting Position
% R2_vec                : same units as R1_vec,                 Ending Position
% ToF_input             : any Time units,                       Time of Flight
% Long                  : true / false,                         Long/Short path, Long being upper branch
% N                     : Integer,                              Number of Revolutions
% mu                    : same units as R1_vec and ToF_input,   Gravitational Parameter of Central Body
% max_iterations        : Integer,                              Number of Iterations

% Outputs               : Units,                                Description
% V1_vec                : same units as R1_vec and ToF_input,   Starting Velocity
% V2_vec                : same units as V1_vec,                 Ending Velocity
% theta                 : deg,                                  Angle between Starting and Ending Position
% semi_major_axis_mid   : same units as R1_vec,                 Semi Major Axis of Transfer Orbit
% eccentricity          : No Units,                             Eccentricity of Transfer Orbit
% DeltaV1               : same units as V_Earth,                Velocity Escaping Earth SOI
% DeltaV2               : same units as V_Mars,                 Velocity Entering Mars SOI

function [DeltaV1, DeltaV2] = Func_Lambert_Mars(R1_vec, R2_vec, ToF_input, Long, N, mu, max_iterations, V_Earth, V_Mars)
%% Setup

r1 = norm(R1_vec);
r2 = norm(R2_vec);
theta = acosd(dot(R1_vec,R2_vec)/(r1*r2));

if theta > 180
    beta_sign = -1;
else
    beta_sign = 1;
end

chord = sqrt(r1*r1 + r2*r2 - 2*r1*r2*cosd(theta)); % Length
semiperimeter = (r1+r2+chord) * 0.5;               % Length

sqrt_mu = sqrt(mu);

semi_major_axis_max = semiperimeter*10000;                     % Length (should be infinity) (10000*semiperimeter)
semi_major_axis_min = semiperimeter/2;             % Length

%% Loop

for i = 1:max_iterations

    % Mid Calcs
    semi_major_axis_mid = (semi_major_axis_min + semi_major_axis_max) * 0.5;                  % Length
    
    if Long == true
        alpha_mid = 2*pi-Func_Alpha(semiperimeter, semi_major_axis_mid);                      % rad
    else
        alpha_mid = Func_Alpha(semiperimeter, semi_major_axis_mid);                           % rad
    end
    
    beta_mid = beta_sign*Func_Beta(semiperimeter, chord, semi_major_axis_mid);                % rad

    ToF_mid = Func_ToF(semi_major_axis_mid, N, alpha_mid, beta_mid, sqrt_mu);                 % Time


    % Max Calcs

    if Long == true
        alpha_max = 2*pi-Func_Alpha(semiperimeter, semi_major_axis_max);                      % rad
    else
        alpha_max = Func_Alpha(semiperimeter, semi_major_axis_max);                           % rad
    end
    
    beta_max = beta_sign*Func_Beta(semiperimeter, chord, semi_major_axis_max);                % rad

    ToF_max = Func_ToF(semi_major_axis_max, N, alpha_max, beta_max, sqrt_mu);                 % Time


    % Min Calcs
    
    if Long == true
        alpha_min = 2*pi-Func_Alpha(semiperimeter, semi_major_axis_min);                      % rad
    else
        alpha_min = Func_Alpha(semiperimeter, semi_major_axis_min);                           % rad
    end
    
    beta_min = beta_sign*Func_Beta(semiperimeter, chord, semi_major_axis_min);                % rad

    ToF_min = Func_ToF(semi_major_axis_min, N, alpha_min, beta_min, sqrt_mu);                 % Time


    % Comparison
    ToF_List_small = [ToF_min, ToF_mid];
    ToF_List_big = [ToF_mid, ToF_max];

    if isbetween(ToF_input,min(ToF_List_small),max(ToF_List_small))
        semi_major_axis_max = semi_major_axis_mid;

    elseif isbetween(ToF_input,min(ToF_List_big),max(ToF_List_big))
        semi_major_axis_min = semi_major_axis_mid;
    else
        % disp('invalid solution (check sma min, max, Long)')
        semi_major_axis_mid = 0;
        break
    end
end

%% Eccentricity Calc

% eccentricity = sqrt(1 - (((4 * (semiperimeter - r1) * (semiperimeter - r2)) / (chord*chord)) * sin((alpha_mid + beta_mid) / 2) * sin((alpha_mid + beta_mid) / 2)));

%% Velocity Calc

A = sqrt(mu/(4*semi_major_axis_mid))*cot(alpha_mid/2);
B = sqrt(mu/(4*semi_major_axis_mid))*cot(beta_mid/2);

U1_vec = R1_vec./r1; 
U2_vec = R2_vec./r2; 
Uc_vec = (R2_vec - R1_vec)./chord;

V1_vec = (B+A)*Uc_vec + (B-A)*U1_vec;
% V1_norm = norm(V1_vec);
V2_vec = (B+A)*Uc_vec - (B-A)*U2_vec;
% V2_norm = norm(V2_vec);

DeltaV1 = norm(V1_vec - V_Earth);
DeltaV2 = norm(V2_vec - V_Mars);

end


function alpha = Func_Alpha(semiperimeter, semi_major_axis)
alpha = 2*asin(sqrt(semiperimeter/(semi_major_axis*2)));
end

function beta = Func_Beta(semiperimeter, chord, semi_major_axis)
beta = 2*asin(sqrt((semiperimeter-chord)/(semi_major_axis*2)));
end

function ToF = Func_ToF(semi_major_axis, N, alpha, beta, sqrt_mu)
ToF = (sqrt(semi_major_axis^3)*(2*N*pi+alpha-beta-(sin(alpha)-sin(beta))))/sqrt_mu;
end

