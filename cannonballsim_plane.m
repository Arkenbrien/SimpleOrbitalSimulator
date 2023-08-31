%% WORKSPACE INIT

close all
clear all
clc

disp('Started the Simulation!')

%% OPTIONS

% time stpes, simulation time (s)
dt          = 1;
t_max       = 10000;

% radius/mass of cannon ball (m)
r_c         = 0.1;
m_c         = 10;

% initial velocity parallel to the surface m/s
v_x(1)      = 1000;
v_z(1)      = 0;

% initial altitude m
H(1)        = 20500;

% drag coefficient of cannonball
Cd          = 0.47;

% Animate??
animate_bool = 1;

%% VAR INIT

% Simulation parameters
m_e         = 5.9722*10^(24);     % mass of earth (kg)
A_c         = pi * r_c^2;         % cross-sectional area (for a sphere) m^2
G           = 6.67408*10^(-11);   % Gravitational Constant m3 kg-1 s-2
r_e         = 6378*10^3;          % radius of earth (m), *assmuing perfect sphere

rho_0       = 1.225;              % Standard air density at sea level (kg/m^3)
L           = 0.0065;             % Temperature lapse rate (K/m)
T_0         = 288.15;             % Standard temperature at sea level (K)
M           = 0.029;              % Molar mass of Earth's air (kg/mol)
R           = 8.314;              % Universal gas constant (J/(mol K))

% Time step and total simulation time
t_array     = 0:dt:t_max;

% initial x coordinates m
x_pos(1)    = 0;

if animate_bool
    figure()
    hold on
end

%% Simulation

%
% Initial Conditions
%

% Calculate gravity
g(1) = (G * m_c * m_e) / (H(1) + r_e)^2 / 10;

% Calculate rho
if H(1) > 20000
    rho = 0;
else
    [~,~,~, rho(1)] = atmosisa(H(1));
end

% Calculate air resistance
F_d_x(1)  = 0.5 * rho(1) * v_x(1)^2 * Cd * A_c;
F_d_z(1)  = 0.5 * rho(1) * v_z(1)^2 * Cd * A_c;

% Calculate Acceleration
a_x(1)    = (-F_d_x * v_x(1)) / m_c;
a_z(1)    = (-m_c * g(1) - F_d_z(1) * v_z(1)) / m_c;

%
% Loop
%

for idx = 2:length(t_array)
    
    % Calculate gravity
    g(idx)      = (G * m_c * m_e) / (H(idx-1) + r_e)^2 / 10;
    
    % Calculate rho
    if H(idx-1) > 20000
        rho(idx) = 0;
    else
        [~,~,~, rho(idx)] = atmosisa(H(idx-1));
    end
    
    % Calculate air resistance
    F_d_x(idx)  = 0.5 * rho(idx) * v_x(idx-1)^2 * Cd * A_c;
    F_d_z(idx)  = 0.5 * rho(idx) * v_z(idx-1)^2 * Cd * A_c;
    
    % Update accelerations
    a_x(idx) = -(F_d_x(idx) / m_c);
    a_z(idx) = -g(idx) + (F_d_z(idx) / m_c);
    
    % Update velocities
    v_x(idx)    = v_x(idx-1) + a_x(idx) * dt;
    v_z(idx)    = v_z(idx-1) + a_z(idx) * dt;
    
    % Update positions
    x_pos(idx)  = x_pos(idx-1) + v_x(idx) * dt;
    H(idx)      = H(idx-1) + v_z(idx) * dt;
    
    % If the cannonball hits the ground, break the loop.
    if H(idx) <= 0
        disp('BOOM')
        break
    end
    
    % DEBUG
%     clc
%     fprintf(' g %f\n rho %f\n Fdx|Fdz %f|%f\n ax|az %f|%f\n vx|vz %f|%f\n x|H %f|%f',...
%             g(idx),rho(idx),F_d_x(idx),F_d_z(idx),a_x(idx),a_z(idx),v_x(idx),v_z(idx),x_pos(idx),H(idx));
        
%     pause(.1)

    if animate_bool
        plot(x_pos(idx), H(idx), 'bo');
        pause(.1)
    end
    
    
    
end

%%
if animate_bool
    scatter(x_pos(idx), H(idx), '*', 'filled', 'MarkerEdgeColor', [1,0,0], 'MarkerFaceColor', [1,0,0], 'SizeData', 1000);
    pause(.5)
    scatter(x_pos(idx), H(idx), '*', 'filled', 'MarkerEdgeColor', [1,.5,0], 'MarkerFaceColor', [1,.5,0], 'SizeData', 500);
    pause(.5)
    scatter(x_pos(idx), H(idx), '*', 'filled', 'MarkerEdgeColor', [1,1,0], 'MarkerFaceColor', [1,1,0], 'SizeData', 200);
end

hold off

%% PLOTZ

figure()
plot(t_array(1:idx),g)
xlabel('time (s)')
ylabel('gravity (m/s^2)')

figure
plot(t_array(1:idx),rho)
xlabel('time (s)')
ylabel('Air Density (kg/m^3)')

figure
plot(x_pos,H)
xlabel('(m)')
ylabel('(m)')

figure
plot(t_array(1:idx),H)
xlabel('time (s)')
ylabel('Altitude (m)')

figure
plot(t_array(1:idx),v_x)
xlabel('time (s)')
ylabel('velocity x (m/s)')

figure
plot(t_array(1:idx),a_z)
xlabel('time (s)')
ylabel('acc z (m/s^2)')

figure
plot(t_array(1:idx),a_x)
xlabel('time (s)')
ylabel('acc x (m/s^2)')

figure
plot(t_array(1:idx),v_z)
xlabel('time (s)')
ylabel('velocity z (m/s)')

figure
plot(t_array(1:idx),F_d_x)
xlabel('time (s)')
ylabel('Drag Force x (N)')

figure
plot(t_array(1:idx),F_d_z)
xlabel('time (s)')
ylabel('Drag Force z (N)')
