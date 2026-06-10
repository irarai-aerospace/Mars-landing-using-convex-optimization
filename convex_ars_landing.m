clc; clear; close all;

%% =======================
%  PS6(f): Convex Mars Landing
%  Requires CVX
%  =======================

%% Parameters from PS6
dt    = 1.0;                    % s
N     = 78;                     % number of intervals
alpha = 0.5086;                 % s/km
g     = [-3.7114e-3; 0; 0];     % km/s^2

Tmin  = 4.97;                   % kg*km/s^2
Tmax  = 13.26;                  % kg*km/s^2

r0    = [1.5; 0; 2.0];          % km
v0    = [-0.075; 0.03; 0.1];    % km/s
m0    = 2000;                   % kg
z0    = log(m0);

gamma_min_deg = 3;              
gamma_min     = deg2rad(gamma_min_deg);
cot_gamma     = 1/tan(gamma_min);

%% alternative
%% Parameters for modified scenario
% Tmin  = 8;
% Tmax  = 20;
% 
% r0    = [2; 4; 6];
% v0    = [-0.078; 0.032; 0.102];
% 
% gamma_min_deg =  3;


%% State and control dimensions
nx = 7;   % x = [r(3); v(3); z]
nu = 4;   % u = [a(3); sigma]

%% Continuous-time matrices from PS6(d.1)
A = [zeros(3,3), eye(3),     zeros(3,1);
     zeros(3,3), zeros(3,3), zeros(3,1);
     zeros(1,3), zeros(1,3), 0];

B = [zeros(3,3), zeros(3,1);
     eye(3),     zeros(3,1);
     zeros(1,3), -alpha];

c = [zeros(3,1);
     g;
     0];

%% Discrete-time matrices (constant for all k)
Ad = [eye(3), dt*eye(3), zeros(3,1);
      zeros(3,3), eye(3), zeros(3,1);
      zeros(1,3), zeros(1,3), 1];

Bd = [0.5*dt^2*eye(3), zeros(3,1);
      dt*eye(3),       zeros(3,1);
      zeros(1,3),      -alpha*dt];

cd = [0.5*dt^2*g;
      dt*g;
      0];

%% Lower bound z_lb(k) = ln(m0 - alpha*Tmax*t_k)
t = 0:dt:N*dt;                  % nodes: k = 0,...,N
zlb = log(m0 - alpha*Tmax*t);   % 1 x (N+1)

%% CVX Optimization
cvx_begin
    cvx_precision best
    variables x(nx, N+1) a(3, N) s(1, N)

    minimize( sum(s) * dt)

    subject to

        % Initial condition
        x(:,1) == [r0; v0; z0];

        % Terminal conditions
        x(1:3, N+1) == [0;0;0];
        x(4:6, N+1) == [0;0;0];
        
        %lateral displacement constraint
        rho_max = 3.5.;   % example in km

        for k = 1:N
          norm( x(2:3,k), 2 ) <= rho_max;
        end

        norm( x(2:3,N+1), 2 ) <= rho_max;

        % Dynamics and path/control constraints
        for k = 1:N

            % control vector u_k = [a_k; sigma_k]
            uk = [a(:,k); s(k)];

            % Discrete dynamics
            x(:,k+1) == Ad*x(:,k) + Bd*uk + cd;

            % z lower bound
            x(7,k) >= zlb(k);

            % Lower thrust bound: Tmin*exp(-z_k) <= sigma_k
            s(k) >= Tmin * exp( -x(7,k) );

            % Upper thrust bound: sigma_k <= Tmax*exp(-zlb_k)*(1 - (z_k - zlb_k))
            s(k) <= Tmax * exp(-zlb(k)) * ( 1 - (x(7,k) - zlb(k)) );

            % SOC thrust/acceleration constraint
            norm( a(:,k), 2 ) <= s(k);

            % Altitude nonnegative: r1 >= 0
            x(1,k) >= 0;

            % Glide-slope:
            % ||[r2; r3]||_2 <= cot(gamma_min)*r1
            norm( x(2:3,k), 2 ) <= cot_gamma * x(1,k);
        end

        % Also enforce z lower bound at final node
        x(7,N+1) >= zlb(N+1);

        % Final node path constraints
        x(1,N+1) >= 0;
        norm( x(2:3,N+1), 2 ) <= cot_gamma * x(1,N+1);

cvx_end

%% Recover useful quantities
r = x(1:3,:);
v = x(4:6,:);
z = x(7,:);

m = exp(z);                     % actual mass history
fuel_consumed = m0 - m(end);    % kg

a_norm = zeros(1,N);
for k = 1:N
    a_norm(k) = norm(a(:,k),2);
end

%% Print results
fprintf('\nCVX status: %s\n', cvx_status);
fprintf('Optimal objective sum(sigma)*dt = %.8f\n', cvx_optval);
fprintf('Final mass = %.8f kg\n', m(end));
fprintf('Consumed fuel = %.8f kg\n', fuel_consumed);

%% Lossless check preview for part (g)
fprintf('Max |sigma - ||a|| | = %.8e\n', max(abs(s - a_norm)));

%% =======================
%  Plots
%  =======================

% 3D trajectory
figure;
plot3(r(1,:), r(2,:), r(3,:), 'LineWidth', 2);
grid on;
xlabel('r_1 (km)');
ylabel('r_2 (km)');
zlabel('r_3 (km)');
title('Mars Landing Trajectory');

% Position history
figure;
plot(t, r(1,:), 'LineWidth', 2); hold on;
plot(t, r(2,:), 'LineWidth', 2);
plot(t, r(3,:), 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Position (km)');
legend('r_1','r_2','r_3');
title('Position History');

% Velocity history
figure;
plot(t, v(1,:), 'LineWidth', 2); hold on;
plot(t, v(2,:), 'LineWidth', 2);
plot(t, v(3,:), 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Velocity (km/s)');
legend('v_1','v_2','v_3');
title('Velocity History');

% Mass history
figure;
plot(t, m, 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Mass (kg)');
title('Mass History');

% sigma and ||a||
figure;
stairs(t(1:end-1), s, 'LineWidth', 2); hold on;
stairs(t(1:end-1), a_norm, '--', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Magnitude');
legend('\sigma_k', '||a_k||_2');
title('Lossless Convexification Check');

% z and z lower bound
figure;
plot(t, z, 'LineWidth', 2); hold on;
plot(t, zlb, '--', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('z = ln(m)');
legend('z','z_{lb}');
title('Log-Mass State and Lower Bound');