%===============================================================================
% ECEN426 Advanced Mechatronic Systems
% Assignment One 2021
%
% This file is skeleton code that implements a Kalman filter to estimate the
% pose of a robot head/neck in the roll dimension.
%
% C.Hollitt
% 20th September 2021
%===============================================================================

%-------------------------------------------------------------------------------
% Form the model.

dt = 10e-3; % 10 ms time interval.
T_max = 20; % Duration of the simulation.
t=[0:1:T_max/dt]'.*dt;

% State vector
% [theta, omega, alpha]

A= [1 dt 0.5*dt^2; 
    0 1 dt; 
    0 0 0];
B= [0 0 1]';
C_1= [0 0 1]; 
C_2= [1 0 0];
D= 0;

% Noise model -------------------

G = [0.5*dt^2 ; dt ; 1];
var_acceleration = 10; % Variance of the acceleration noise on the head.

Q = G * var_acceleration * (G');

R_1 = 1;   % Variance in the acceleration measurement (sensor 1).
R_2 = 0.5; % Variance in the position measurement (sensor 2).

%-------------------------------------------------------------------------------
% Initial  Conditions.

% Actual position that the neck begins the simulation.
x0 = [-2, 0, 10]'; % Don't change this - it is the true initial state that you 
                  % "don't know".

x(:,1) = x0; % Set the initial state for the simulation.

% Our initial guess at the neck position. a)
x_post(:,1) = [0,0,0]';     % Put your initial estimate here. 
P_post = diag([5,5,5]);     % Put your initial estimate here. 

% %-----------------------------------------------------------------------------
L_acc = zeros(length(t),3);
L_acc(1,:) = [P_post * C_1' * inv( C_1 * P_post *C_1' + R_1)].';  % initialise the accumulated L with initial guess
% % Run simulation b)

I=eye(3);

for k=2:length(t)
  % x(k) = Ax(k-1) + Bu(k) + w(k)
  w = G*sqrt(var_acceleration)*randn(1,1);

  x(:,k) = A*x(:,k-1) + w;   % u=0, so omit the B*u(k) term.
  
  x_prior(:,k) = A*x_post(:,k-1);
  P_prior = A*P_post*A' + Q;
  
  % Update the state using Sensor One at every step. 
  % If you were to change S1_update_interval to 3 for example, then this block
  % of code would only run every third time through the outer loop.

  S1_update_interval = 1; 
  if mod(k, S1_update_interval)==0 
    v = sqrt(R_1)*randn(1,1); % Find noise to add to the measurement.
    y(:,k) = C_1*x(:,k) + v;         % Measure the state.

    L = P_prior * C_1' * inv( C_1 * P_prior *C_1' + R_1);
    L_acc(k,:) = L;
    x_post(:,k) = x_prior(:,k) + L*(y(:,k)-C_1*x_prior(:,k));
    P_post = (I-L*C_1)*P_prior;
    P_post = 0.5*(P_post + P_post'); % Make sure that the covariance remains symmetric.
  end
  
end

% Plot results -----------------------------------------------------------------
figure(1)
for plot_num = 1:3
  subplot(3,1,plot_num)
  stairs(t,[x(plot_num,:);x_post(plot_num,:)]')
  ylabel(['x_' num2str(plot_num)])
  legend("x","x posterior")
end
xlabel('Time (s)')

% B
figure(2)
plot(t,L_acc)
ylabel("L")
xlabel("Time (s)")
legend("x1","x2","x3")
title("Startup Transients")

%C Steady State Kalman Filter

%L is constant and set to the last value of L in the non-steady state run
%through

L_ss = L;
x_ss(:,1) = x0; % Set the initial state for the simulation.
x_post_ss(:,1) = [0,0,0]';     % Put your initial estimate here. 

for k=2:length(t)
  % x(k) = Ax(k-1) + Bu(k) + w(k)
  w = G*sqrt(var_acceleration)*randn(1,1);

  x_ss(:,k) = A*x_ss(:,k-1) + w;   % u=0, so omit the B*u(k) term.
  
  x_prior_ss(:,k) = A*x_post_ss(:,k-1);
  
  % Update the state using Sensor One at every step. 
  % If you were to change S1_update_interval to 3 for example, then this block
  % of code would only run every third time through the outer loop.

  S1_update_interval = 1; 
  if mod(k, S1_update_interval)==0 
    v = sqrt(R_1)*randn(1,1); % Find noise to add to the measurement.
    y(:,k) = C_1*x_ss(:,k) + v;         % Measure the state.

    x_post_ss(:,k) = x_prior_ss(:,k) + L_ss*(y(:,k)-C_1*x_prior_ss(:,k));
  end
  
end

% Plot results -----------------------------------------------------------------
figure(3)
for plot_num = 1:3
  subplot(3,1,plot_num)
  stairs(t,[x_ss(plot_num,:);x_post_ss(plot_num,:)]')
  ylabel(['x_' num2str(plot_num)])
  legend("x","x posterior")
end
xlabel('Time (s)')
%Is the same because the non-steady-state filter quickly settles on L

% D Add a second sensor

x_sen(:,1) = x0; % Set the initial state for the simulation.
% Our initial guess at the neck position. a)
x_post_sen(:,1) = [0,0,0]';     % Put your initial estimate here. 
P_post_sen = diag([5,5,5]);     % Put your initial estimate here. 

for k=2:length(t)
  % x(k) = Ax(k-1) + Bu(k) + w(k)
  w = G*sqrt(var_acceleration)*randn(1,1);

  x_sen(:,k) = A*x_sen(:,k-1) + w;   % u=0, so omit the B*u(k) term.
  
  x_prior_sen(:,k) = A*x_post_sen(:,k-1);
  P_prior_sen = A*P_post_sen*A' + Q;
  
  % Update the state using Sensor One at every step. 
  % If you were to change S1_update_interval to 3 for example, then this block
  % of code would only run every third time through the outer loop.

  S1_update_interval = 1; 
  S2_update_interval = 30;
  if mod(k, S1_update_interval)==0 
    v = sqrt(R_1)*randn(1,1); % Find noise to add to the measurement.
    if mod(k, S2_update_interval) == 0
        v2 = sqrt(R_2)*randn(1,1);
        
        y(:,k) = C_2*x_sen(:,k) + v2;         % Measure the state.
        L_sen = P_prior_sen * C_2' * inv( C_2 * P_prior_sen *C_2' + R_2);
        L_acc_sen(k,:) = L_sen;
        x_post_sen(:,k) = x_prior_sen(:,k) + L_sen*(y(:,k)-C_2*x_prior_sen(:,k));
        P_post_sen = (I-L_sen*C_2)*P_prior_sen;
        P_post_sen = 0.5*(P_post_sen + P_post_sen'); % Make sure that the covariance remains symmetric.

    else
        y(:,k) = C_1*x_sen(:,k) + v;         % Measure the state.

        L_sen = P_prior_sen * C_1' * inv( C_1 * P_prior_sen *C_1' + R_1);
        L_acc_sen(k,:) = L_sen;
        x_post_sen(:,k) = x_prior_sen(:,k) + L_sen*(y(:,k)-C_1*x_prior_sen(:,k));
        P_post_sen = (I-L_sen*C_1)*P_prior_sen;
        P_post_sen = 0.5*(P_post_sen + P_post_sen'); % Make sure that the covariance remains symmetric.
    end
  end
  
  
  
end
figure(4)
for plot_num = 1:3
  subplot(3,1,plot_num)
  stairs(t,[x_sen(plot_num,:);x_post_sen(plot_num,:)]')
  ylabel(['x_' num2str(plot_num)])
  legend("x","x posterior")
end
xlabel('Time (s)')

% END ==========================================================================