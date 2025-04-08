function [x, y] = covariance_ellipse(P,ux,uy,CEP)
% [x y] = covariance_ellipse(P,ux,uy,CEP): 	
%			This function plots the covariance ellipse
%			of a distribution with a 2x2 covariance
%			matrix P and (optional) means of ux and uy.
%			(Default: ux=uy=0)
%
%			You can specify CEP as either 50 or 90, meaning
%			that an ellipse will be created such that 50% or 90%
%			of all data points in the given distribution will
%			lie within the ellipse.  (Default: CEP=90)
%
%			Function will plot ellipse unless [x y] is specified.


% Perform Input Error Checking
	if nargin < 1
   	error('You must input a covariance matrix')
	end
	[m n]=size(P);
	if m~=2 | n~=2
   	error('Inputted covariance matrix must be 2x2')
	end
% 	if abs(P(1,2)- P(2,1)) > .01*abs(P(1,2))
%    	error('Inputted covariance matrix must be symmetrical')
% 	end
	if nargin == 1
   	ux = 0;
	   uy = 0;
	end
	if nargin >= 2 & length(ux) ~= 1
   	error('Mean of x should be a single number')
	end
	if nargin >= 3 & length(uy) ~= 1
   	error('Mean of y should be a single number')
   end
   if nargin >= 4 & length(CEP) ~= 1
      error('Error Probability (CEP) should be a single number')
   end
   if nargin < 4
      CEP = 90;
   end 
% End Input Error Checking


% Covariance Matrix is of form:
%		[ sigX^2				rho*sigX*sigY]
%		[ rho*sigX*sigY	sigY^2       ]

% Find X and Y standard deviations.  Find rho.
sigX = sqrt(P(1,1));
sigY = sqrt(P(2,2));
rho  = P(2,1)/sigX/sigY;

% When k=1.39, the probability that a point (x,y) is within the ellipse is 50%
% When k=4.61, the probability that a point (x,y) is within the ellipse is 90%
if CEP == 90
   k = 4.61; 	
elseif CEP == 50
   k = 1.39;
else
   disp('Warning: Invalid CEP.  CEP set to 90')
   CEP = 90;
   k = 4.61;
end

% number of points in ellipse (must be even)
npoints=500; 

% horizontal points of ellipse (with mean of 0)
xbound = sqrt(k*sigX^2);
x0 = [linspace(-xbound,xbound,npoints/2) linspace(xbound,-xbound,npoints/2)];

% vertical points of ellipse (with mean of 0)
idx = 1:npoints/2;  % relevant indices
y0(idx) = rho*sigY/sigX*x0(idx) - sigY*sqrt( (x0(idx).^2/sigX^2-k)*(rho^2-1) );
idx = (1:npoints/2)+npoints/2;
y0(idx) = rho*sigY/sigX*x0(idx) + sigY*sqrt( (x0(idx).^2/sigX^2-k)*(rho^2-1) );

% Center ellipse at (ux,uy)
x = x0 + ux;
y = real(y0 + uy);	% y may have SMALL imaginary part.  Cut it.

% Plot ellipse if there are no output arguments
if nargout == 0
	plot(x,y,ux,uy,'rx')
	plot_width = max([abs(ux)+sqrt(5)*sigX abs(uy)+sqrt(5)*sigY]);
	axis([-plot_width plot_width -plot_width plot_width])
	axis square
   grid
   title('Covariance Ellipse')
   xlabel('X units')
   ylabel('Y units')
end

