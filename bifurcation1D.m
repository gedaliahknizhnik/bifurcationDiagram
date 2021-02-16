function [bfData] = bifurcation1D(dynSys, xLim, rLim)
% BIFURCATION1D This function creates a bifurcation diagram for a 1D
% system given its dynamics in the form xdot = f(x,r).
%
% INPUTS:
%   dynSys - [required] function handle to a function of the form
%               [xdot, fprime] = dynSys(x,r). 
%            Function must be vectorized to return values for multiple
%            values of x.
%   xLim   - [optional] 1x1 max value of x. Will consider -xLim<x<xLim
%   rLim   - [optional] 1x1 max value of r. Will consider -rLim<r<rLim
%
% OUTPUTS:
%   bfData - nx3 where each row is a fixed point of the form
%               [r, xc, stability] where xc is the critical point and
%               stability is -1 for stable, 0 for half-stable, and 1 for
%               unstable.


if nargin < 2
    xLim = pi/2;
end

if nargin < 3
    rLim = 4;
end
    xs = [-xLim:0.01:xLim]';

bfData = [];

for rr = -rLim:0.05:rLim

    [dxs, fPrimes] = dynSys(xs,rr);

    [zs,ss] = zeroFinder(dxs, xs, fPrimes);

    bfData = [bfData; rr*ones(size(zs,1),1), zs, ss];
end
close all

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


function [zs, ss] = zeroFinder(y, x, f)
    
    % Find spots actually equal to zero
    zeroInds = find(y == 0);
    zs       = x(zeroInds);
    ss       = sign(f(zeroInds));
    
    % Find zero crossings
    zci = find(y(:).*circshift(y(:), [-1 0]) <= 0);  
    
    
    % Remove zero crossings taken care of by zeroInds
    for ii=1:size(zeroInds,1)
        zci(abs(zci(:,1) - zeroInds(ii)) <= 1) = [];
    end
    
    if numel(zci) > 0
        % Remove last index if needed
        if zci(end) == size(y,1)
            zci(end) = [];
        end
        
        for ii=1:size(zci,1)
            zs = [zs; interp1(y(zci(ii):zci(ii)+1),x(zci(ii):zci(ii)+1),0)];
            ss = [ss; sign(interp1(y(zci(ii):zci(ii)+1),f(zci(ii):zci(ii)+1),0))];
        end
    end
    zs;
end