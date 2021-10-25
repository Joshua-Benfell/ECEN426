clear
clf
clc

% Set up the environment (the landmark locations).
x_beacons = [20 20  40 20  60 20  80 20  20 -10  40 -10  60 -10  80 -10]';
%x_beacons = [20 20  40 20]';
P_beacons = 0.001^2*eye(length(x_beacons));
%P_beacons = 10^2*eye(length(x_beacons));

x_robot = [0 0 0]';     % The robot's true starting position.

% Robot process noise.
sigma_forward = 1; % Standard deviation of noise in forward direction [m].
                   %  This should really be related to u1, because the
                   %  noise would depend on the robot's motion, and
                   %  certainly go to zero if the robot were static.
sigma_bearing = 2; % Standard deviation of turning noise [deg].
Q = diag([sigma_forward^2 sigma_bearing^2]);

% Measurement noise.
sigma_y_range = 1;      % Std dev of range measurements [m]
sigma_y_bearing = 2;    % Std dev of angle bearing measurements [deg]
R = diag([sigma_y_range^2 sigma_y_bearing^2]);

% Define the full state vector.
global xest P
xest = [nan; nan; nan; x_beacons];
P = blkdiag(0*eye(3), P_beacons);
initialise_localiser();

%plot_estimates(xest, P);
xlim([-20 120])
ylim([-40 40])
axis equal

for t = 1:30
    % Move the real robot.
    
    [u1, u2] = controller_command(t, xest(1:3));
    
    w = mvnrnd([0 0], Q); % Process noise.
    %w = [0 0];
    
    x_robot(1) = x_robot(1) + (u1+w(1))*cos(deg2rad(x_robot(3)));
    x_robot(2) = x_robot(2) + (u1+w(1))*sin(deg2rad(x_robot(3)));
    x_robot(3) = x_robot(3) + u2 + w(2);
    hold on
    plot(x_robot(1), x_robot(2), 'm+')
    hold off
    
    predict([u1 u2]);
    
    % Generate the measurements.
    
    v = mvnrnd([0 0], R); % Measurement noise.
    %v = [0 0];
    
    for beacon = 1:2:length(x_beacons)
      y(beacon,:) = sqrt((x_robot(1) - x_beacons(beacon))^2 + ...
                        (x_robot(2) - x_beacons(beacon+1))^2 ) + v(1);
      y(beacon+1,:) = rad2deg(atan2((x_robot(2) - x_beacons(beacon+1)),...
                          (x_robot(1) - x_beacons(beacon)))) ...
                          - x_robot(3) + v(2);
    end
    
    hold on
    
    % Take a measurement only every meas_interval-th time step.
    meas_interval = 1;
    if mod(t, meas_interval) == 0
        correct(y);
    end

    pause(0.25)    
end

%==========================================================================
% Contoller Code
%==========================================================================

function [u1, u2] = controller_command(t, x)
% Execute a preset plan for the robot motion. 
% The robot's state (x) is not needed for this particular plan.
    u1_plan = [10 10 10 10 10 10 10 ...
               5 5 5 5 5 5 5 ...
               5 5 5 5 5 5 ...
               5 5 5 5 5 5 ...
               10 10 10 10 10 10];
    u2_plan = [0 0 0 0 0 0 0 ...
               15 15 15 15 15 15 15 ...
               0 0 0 0 0 0 ...
               15 15 15 15 15 15 ...
               0 0 0 0 0 0];
    
    if t<length(u1_plan)
        u1 = u1_plan(t);
        u2 = u2_plan(t);
    else
        u1 = 0;
        u2 = 0;
    end
end

%==========================================================================
% Localisation code
%==========================================================================


%--------------------------------------------------------------------------
% The Prediction Step.
%--------------------------------------------------------------------------

function [] = initialise_localiser()
    global xest P

    xest(1:3) = [0 0 0]';
    P(1,1) = 1^2;   % Initial variance of x-position estimate.
    P(2,2) = 1^2;   % Initial variance of y-position estimate.
    P(3,3) = 5^2;   % Initial variance of heading estimate.
end

function [] = predict(u) 
    global xest P
    
    sigma_forward = 1;
    sigma_bearing = 2;
    Q = diag([sigma_forward^2 sigma_bearing^2]);

    % Form estimate of robot's position.
    xest(1) = xest(1) + u(1)*cos(deg2rad(xest(3)));
    xest(2) = xest(2) + u(1)*sin(deg2rad(xest(3)));
    xest(3) = xest(3) + u(2);
    
    % Update the convariance of our estimate of the robot's position.
    A = eye(3);
    G = [cos(deg2rad(xest(3))) 0;
         sin(deg2rad(xest(3))) 0
         0 1];
    P(1:3, 1:3) = A*P(1:3, 1:3)*A' + G*Q*G';
   
    hold on
    plot_estimates(xest, P, 2, 'g', 'r');
end

function [] = correct(y)
    global xest P
    I = eye(3 + length(y));
    
    % Create the noise model. Assume all measurements have the same noise.
    sigma_y_range = 2;      % Std dev of range measurements [m]
    sigma_y_bearing = 2;    % Std dev of angle bearing measurements [deg]
    R = diag(repmat([sigma_y_range^2 sigma_y_bearing^2], [1,length(y)/2]));
    
    % Infer the appropriate C matrix from the size of y.
    for beacon = 1:2:length(y/2)
        C(beacon,:) = [(xest(1) - xest(3 + beacon)) / y(beacon) ...
                       (xest(2) - xest(4 + beacon)) / y(beacon) ...
                       0 zeros(1, length(P) - 3)];
        C(beacon+1,:) = [-(xest(2) - xest(4 + beacon)) / y(beacon)^2 ...
                          (xest(1) - xest(3 + beacon)) / y(beacon)^2 ...
                        -1 zeros(1, length(P) - 3)];
        C(beacon, beacon+3) = -C(beacon,1);
        C(beacon, beacon+4) = -C(beacon,2);
        C(beacon+1, beacon+3) = -C(beacon+1,1);;
        C(beacon+1, beacon+4) = -C(beacon+1,1);;
    end
    
    % Find the optimal Kalman gain.
    L = P * C' * inv(C*P*C' + R);
    
    % Find the expected outputs given the model.
    for beacon = 1:2:length(y/2)
      yest(beacon,:) = sqrt((xest(1) - xest(beacon+3))^2 + ...
                        (xest(2) - xest(beacon+4))^2 );
      yest(beacon+1,:) = rad2deg(atan2((xest(2) - xest(beacon+4)),...
                          (xest(1) - xest(beacon+3)))) - xest(3);
    end
    
    % Update the estimate.
    xest = xest + L*(y-yest);
    P = (I-L*C)*P;
    
    hold on
    plot_estimates(xest, P, 2, 'b', 'r');
end