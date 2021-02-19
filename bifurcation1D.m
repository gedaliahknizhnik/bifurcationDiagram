function [bfData] = bifurcation1D(dynSys, xLim, rLim)
% BIFURCATION1D This function creates a bifurcation diagram for a 1D
%   system given its dynamics in the form xdot = f(x,r).
%
%   The function will handle inputs that are vectorized or not. For best
%   performance, vectorize in BOTH x and r to allow a meshgrid operation to
%   succeed. Ok performance is achieved with vectorization in x only. Worst
%   performance occurs if not vectorized at all.
%
% INPUTS:
%   dynSys - [required] function handle to a function of the form
%               [xdot, fprime] = dynSys(x,r). 
%   xLim   - [optional] 1x1 max value of x. Will consider -xLim<x<xLim
%   rLim   - [optional] 1x1 max value of r. Will consider -rLim<r<rLim
%
% OUTPUTS:
%   bfData - nx3 where each row is a fixed point of the form
%               [r, xc, stability] where xc is the critical point and
%               stability is -1 for stable, 0 for half-stable, and 1 for
%               unstable.


if nargin < 2
    xLim = 3.5;%pi/2;
end

if nargin < 3
    rLim = 4;
end

step = 0.01;

xs = [-xLim:step:xLim]';
rs = [-rLim:step:rLim]';

bfData = [];

% Attempt to solve using meshgrid for optimal time
% Function needs to be vectorized in both x and r
try
    [Xs,Rs]        = meshgrid(xs,rs); % Assemble the input data
    [dXs, fPrimes] = dynSys(Xs,Rs);   % Evaluate the dynamic system
    [zs,ss,rs]     = zeroFinder(dXs, Xs, fPrimes, Rs);    % Find the zeros
    bfData         = [rs,zs,ss];      % Assemble the final data

% If function isn't vectorized in both x and r, try 
% to use x as a vector and r as a scalar
catch
    try
        % For each r
        for rr = -rLim:step:rLim
            [dxs, fPrimes] = dynSys(xs,rr); % Evaluate the system
            [zs,ss] = zeroFinder(dxs, xs, fPrimes, rr); % Find the zeros
            bfData = [bfData; rr*ones(size(zs,1),1), zs, ss]; % Assemble the final data
        end
        
   % If the function isn't vectorized in either x or r, go through each one
    catch
        
        % For each r
        for rr = -rLim:step:rLim
            
            dxs     = 0*xs;
            fPrimes = dxs;
            
            % For each x, assemble the dxs and fPrimes
            for ii=1:numel(xs)
                [dx, fPrime] = dynSys(xs(ii), rr); % Evaluate the system
                dxs(ii)      = dx;
                fPrimes(ii)  = fPrime;
            end
            
            [zs,ss] = zeroFinder(dxs, xs, fPrimes, rr);       % Find the zeros
            bfData = [bfData; rr*ones(size(zs,1),1), zs, ss]; % Assemble the final data
        end
    end
end
    
%% Plot the bifurcation diagram
figure

iStable  = find(bfData(:,3) < 0);
iUstable = find(bfData(:,3) > 0);
iMarg    = find(bfData(:,3) == 0);

subplot(1,2,1), hold on, grid on
plot(bfData(iStable,1),  bfData(iStable,2) ,'b*','Linewidth',1.2)
plot(bfData(iUstable,1), bfData(iUstable,2),'r*','Linewidth',1.2)
plot(bfData(iMarg,1),    bfData(iMarg,2)    ,'g*','Linewidth',1.2)
xlabel('$r$','Interpreter','Latex')
ylabel('$x$','Interpreter','Latex')
xlim([-rLim,rLim]), axis square

subplot(1,2,2), hold on, grid on
plot(bfData(iStable,2),  bfData(iStable,1) ,'b*','Linewidth',1.2)
plot(bfData(iUstable,2), bfData(iUstable,1),'r*','Linewidth',1.2)
plot(bfData(iMarg,2),    bfData(iMarg,1)    ,'g*','Linewidth',1.2)
xlabel('$x$','Interpreter','Latex')
ylabel('$r$','Interpreter','Latex')
ylim([-rLim,rLim]), axis square

end