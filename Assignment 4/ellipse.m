function [x_out,y_out] = ellipse(A, b, n_points, dimensions, style)
% Plot an ellipse defined by xA'x' + b = 0
% A is the covariance matrix.
% If A is larger than R^{2*2} then "dimensions" can be used to select the
%   elements that are to be plotted.
%   For example, dimensions = [1,4] will plot an ellipse of x_1 vs x_4
% style is the linestyle used for the plot.
%--------------------------------------------------------------------------

%==========================================================================
% $Rev: 14 $
% $Date: 2014-09-17 (Wed, 17 Sep 2014) $
% $Author: hollitch $
%==========================================================================

if nargin < 2
  b=[0,0]';
end

if nargin < 3
  n_points = 100;
end

if nargin < 5
  style= 'b';
end

if nargin > 3
  i1 = dimensions(1);
  i2 = dimensions(2);
  A_temp(1,1)=A(i1,i1);
  A_temp(1,2)=A(i1,i2);
  A_temp(2,2)=A(i2,i2);
  A_temp(2,1)=A(i2,i1);
  A=A_temp;
  
  if length(b) > 2 % If the length is 2 then the centre is already in the selected dimensions.
    b_temp(1) = b(i1);
    b_temp(2) = b(i2);
    b= b_temp;
  end
end

theta = [0:2*pi/n_points:2*pi];
x_1 = cos(theta);
x_2 = sin(theta);
x=[x_1 ; x_2];

for ii = 1:length(x);
  radius = [x(1,ii);x(2,ii)]' * inv(A) * [x(1,ii);x(2,ii)];
  y(:,ii) = 1./sqrt(radius) * x(:,ii);
end

x_out = y(1,:) + b(1);
y_out = y(2,:) + b(2);

if nargout == 0 
  plot(x_out, y_out, style)
end