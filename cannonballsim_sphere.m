%% WORKSPACE INIT

close all
clear all
clc

disp('Started the Simulation!')

%% OPTIONS

% time stpes
dt          = .5;

% Num time-steps before exiting due to the cannonball not being able to
% reach the ground
idx_max     = 1000;

% radius/mass of cannon ball (m)
r_c         = 0.1;
m_c         = 10;

% initial velocity parallel to the surface m/s
v_x(1)      = 3000;
v_z(1)      = 0;

% initial altitude m
H(1)        = 30000;

% drag coefficient of cannonball
Cd          = 0.47;

% Animate?
ani_bool    = 1;

% wait until 

%% VAR INIT

idx         = 2;                    % indexes start at 2 >:)

% Simulation parameters
m_e         = 5.9722*10^(24);       % mass of earth (kg)
A_c         = pi * r_c^2;           % cross-sectional area (for a sphere) m^2
G           = 6.67408*10^(-11);     % Gravitational Constant (m^3/kg*s)
r_e         = 6378*10^3;            % radius of earth (m), *assmuing perfect sphere

rho_0       = 1.225;                % Standard air density at sea level (kg/m^3)
L           = 0.0065;               % Temperature lapse rate (K/m)
T_0         = 288.15;               % Standard temperature at sea level (K)
M           = 0.029;                % Molar mass of Earth's air (kg/mol)
R           = 8.314;                % Universal gas constant (J/(mol K))

% initial x/y coordinates m
x_pos(1)    = 0;
z_pos(1)    = H(1) + r_e;

% Animation setup
if ani_bool
    boom = false; % if it doesn't hit the ground it'll fail. :(
    figure();
    hold on
    viscircles([0 0], r_e, 'Color', 'k') % surface of sphere
    viscircles([0 0], r_e+20000, 'Color', 'c') % atmospheric barrier
    axis equal
end

%% Simulation

%
% Initial Conditions
%

% Calculate gravity
g(1)        = (G * m_c * m_e) / (H(1) + r_e)^2 / 10;

% Calculate rho
if H(1) > 20000
    rho         = 0;
elseif H(1) < 20000
    [~,~,~, rho(1)] = atmosisa(H(1));
end

% Calculate air resistance
F_d_x(1)    = 0.5 * rho(1) * v_x(1)^2 * Cd * A_c;
F_d_z(1)    = 0.5 * rho(1) * v_z(1)^2 * Cd * A_c;

% Calculate Acceleration
a_x(1)      = (-F_d_x * v_x(1)) / m_c;
a_z(1)      = (-m_c * g(1) - F_d_z(1) * v_z(1)) / m_c;

%
% Loop
%

while true
    
    % Calculate gravity
    g(idx)          = (G * m_c * m_e) / (H(idx-1) + r_e)^2 / 10;
    p_theta(idx)    = atan2(z_pos(idx-1), x_pos(idx-1)); 
    g_x(idx)        = g(idx)*cos(p_theta(idx));
    g_z(idx)        = g(idx)*sin(p_theta(idx));
    

    
    % Calculate rho
    if H(idx-1) > 20000
        rho(idx)    = 0;
    else
        [~,~,~, rho(idx)] = atmosisa(H(idx-1));
    end
    
    % Calculate air resistance
    v(idx)          = sqrt(v_x(idx-1)^2 + v_z(idx-1)^2);
    v_theta(idx)    = atan2(v_z(idx-1), v_x(idx-1));
    
    F_d(idx)        = 0.5 * rho(idx) * v(idx)^2 * Cd * A_c;
    F_d_x(idx)      = F_d(idx)*cos(v_theta(idx))*-1;
    F_d_z(idx)      = F_d(idx)*sin(v_theta(idx))*-1;
    
    % stupid test - leave it???
    if x_pos >=0
        g_x(idx) = g_x(idx) * -1;
    end
    if z_pos >= 0
        g_z(idx) = g_z(idx) * -1;
    end
    
    % Update Accelerations
    a_x(idx) = g_x(idx) + (F_d_x(idx) / m_c);
    a_z(idx) = g_z(idx) + (F_d_z(idx) / m_c);
    
    % Update velocities
    v_x(idx)    = v_x(idx-1) + a_x(idx) * dt;
    v_z(idx)    = v_z(idx-1) + a_z(idx) * dt;
    
    % Update Positions
    x_pos(idx)  = x_pos(idx-1) + v_x(idx) * dt;
    z_pos(idx)  = z_pos(idx-1) + v_z(idx) * dt;
    H(idx)      = sqrt(x_pos(idx)^2 + z_pos(idx)^2) - r_e;
    
    % If the cannonball hits the ground, break the loop.
    if H(idx) <= 0
        boom = true;
        disp('BOOM')
        break
    end
    
    % DEBUG
%     clc
%     fprintf(' g %f\n rho %f\n Fdx|Fdz %f|%f\n ax|az %f|%f\n vx|vz %f|%f\n x|H %f|%f',...
%             g(idx),rho(idx),F_d_x(idx),F_d_z(idx),a_x(idx),a_z(idx),v_x(idx),v_z(idx),x_pos(idx),H(idx));
%         
%     pause(.1)

    if ani_bool
        plot(x_pos(idx), z_pos(idx), 'bo');
        pause(.1)
    end
    
    if idx > idx_max
        disp('Too fast! Exiting orbit!')
        break
    end
    
    % increase the counter
    idx = idx + 1;
    
end

if ani_bool && boom
    scatter(x_pos(idx), z_pos(idx), '*', 'filled', 'MarkerEdgeColor', [1,0,0], 'MarkerFaceColor', [1,0,0], 'SizeData', 1000);
    pause(.5)
    scatter(x_pos(idx), z_pos(idx), '*', 'filled', 'MarkerEdgeColor', [1,.5,0], 'MarkerFaceColor', [1,.5,0], 'SizeData', 500);
    pause(.5)
    scatter(x_pos(idx), z_pos(idx), '*', 'filled', 'MarkerEdgeColor', [1,1,0], 'MarkerFaceColor', [1,1,0], 'SizeData', 200);
end

hold off

%% PLOTZ

% Time Array
t_array = linspace(0,dt,length(g));

if ~ani_bool
    figure()
    hold on
    viscircles([0 0], r_e, 'Color', 'k') % surface of sphere
    viscircles([0 0], r_e+20000, 'Color', 'c') % atmospheric barrier
    scatter(x_pos, z_pos)
    scatter(x_pos(idx), z_pos(idx), '*', 'filled', 'MarkerEdgeColor', [1,0,0], 'MarkerFaceColor', [1,0,0], 'SizeData', 1000);
    scatter(x_pos(idx), z_pos(idx), '*', 'filled', 'MarkerEdgeColor', [1,.5,0], 'MarkerFaceColor', [1,.5,0], 'SizeData', 500);
    scatter(x_pos(idx), z_pos(idx), '*', 'filled', 'MarkerEdgeColor', [1,1,0], 'MarkerFaceColor', [1,1,0], 'SizeData', 200);
    axis equal
end

figure(); hold on;
plot(t_array,g)
plot(t_array,g_x)
plot(t_array,g_z)
xlabel('time (s)')
ylabel('gravity (m/s^2)')

figure
plot(t_array,rho)
xlabel('time (s)')
ylabel('Air Density (kg/m^3)')

figure
plot(t_array,H)
xlabel('time (s)')
ylabel('Altitude (m)')

figure; hold on
plot(t_array,a_z)
plot(t_array,a_x)
xlabel('time (s)')
ylabel('acc z (m/s^2)')
legend('Z', 'X')
hold off

figure; hold on
plot(t_array,v_z)
plot(t_array,v_x)
xlabel('time (s)')
ylabel('velocity (m/s)')
legend('Z', 'X')
hold off

figure
plot(t_array,F_d_x)
xlabel('time (s)')
ylabel('Drag Force x (N)')

figure
plot(t_array,F_d_z)
xlabel('time (s)')
ylabel('Drag Force z (N)')

figure; hold on
plot(t_array, v_theta)
plot(t_array, p_theta)
xlabel('time (s)')
ylabel('Angle (rad)')
legend({'velocity', 'position'})



