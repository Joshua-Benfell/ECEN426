function [] = plot_estimates(xest, P, show_heading, robot_colour, landmark_colour)
% Plot position and variance in estiamtes of locations.

% The first three elements of xest are the estimated x,y and heading of the
% robot.
% Subsequent pairs of state variables are the location of a landmark.
%   eg x4 is the x locatio of landmark 1 and x5, is its y location.
% if show_heading == 0, then only the error ellipse for the robot's position
%                      is dispplayed.
% if show_heading > 0, then a line showing the direction of estimated robot
%                      travel is shown. Its magnitude is meaningless.
% if show_heading == 2 the an extra two lines are shown, showing the extent
%                      of a one standard deviation error in the estimate.
% robot_colour and landmark_colour use the normal matlab convention for
% defining colours, such as 'b' or [0 0 1].


    if nargin < 3
        show_heading = 2;
    end

    if nargin < 5 
        robot_colour = 'b';
    end
    ellipse(P, xest, 64, [1 2], robot_colour);
    plot(xest(1), xest(2), [robot_colour '+']);
    
    % Plot a wedge showing the estimated heading.
    hold on
    
    if show_heading
        scale = sqrt(max(P(1,1), P(2,2)));
        heading = deg2rad(xest(3));
        line([xest(1) xest(1) + scale * cos(heading)], ...
             [xest(2) xest(2) + scale * sin(heading)], 'Color', robot_colour);
    end
        
    if show_heading == 2
        ang = deg2rad(sqrt(P(3,3)));
        line([xest(1) xest(1) + scale * cos(heading+ang)], ...
             [xest(2) xest(2) + scale * sin(heading+ang)], 'Color', robot_colour);
        line([xest(1) xest(1) + scale * cos(heading-ang)], ...
             [xest(2) xest(2) + scale * sin(heading-ang)], 'Color', robot_colour);
    end
     
    % Plot the positions of landmarks.
    if length(xest) > 3
        if nargin < 5 
            landmark_colour = 'r';
        end
        for this_index = 4:2:length(xest)
            ellipse(P, xest, 64, [this_index, this_index+1], landmark_colour);
            plot(xest(this_index), xest(this_index+1),...
                 [landmark_colour '+']);
        end
    end
    hold off
end

