function [zs, ss, rs] = zeroFinder(y, x, f, r)
% ZEROFINDER This function finds zeros in an input function of the form
%   y = g(x) and returns their stability as a function of g'(x). This is
%   done for a bifurcation diagram, so involves the parameter r.
%
%   This function operates in two modes - meshgrid and vector. The inputs
%   are passed as either meshgrids of comparable dimension or as column
%   vectors.
%
% INPUTS:
%   y - [MxN] or {Mx1] the value of the function y = g(x)
%   x - [MxN] or [Mx1] the inputs to the function y = g(x)
%   f - [MxN] or [Mx1] the value of the derivative of the function g'(x)
%   r - [MxN] or [1x1] the parameter r
%
% OUTPUTS:
%   zs - [Kx1] vector of values of x at which the function is zero.
%   ss - [Kx1] vector of stabilities of corresponding zeros (-1 for stable,
%               0 for marginal, and 1 for unstable)
%   rs - [Kx1] vector of r parameters corresponding to the zeros.

    % Find spots actually equal to zero
    zeroInds = find(y == 0);
    zs       = x(zeroInds);
    ss       = sign(f(zeroInds));
    
    if numel(r) > 1
        rs = r(zeroInds);
    else
        rs       = r + 0*zs;
    end
    
    % Find zero crossings
    zci = find(y(:).*circshift(y(:), [-1 0]) <= 0);  
    
    % Remove zero crossings taken care of by zeroInds
    for ii=1:size(zeroInds,1)
        zci(abs(zci(:,1) - zeroInds(ii)) <= 1) = [];
    end
    
    % Remove entries that correspond to ends of columns;
    zci(mod(zci,size(x,1)) == 0) = [];        
    
    % If there are any elements left
    if numel(zci) > 0
        % Interpolate to find the zero value and add it and its stability        
        for ii=1:size(zci,1)
            zs = [zs; interp1(y(zci(ii):zci(ii)+1),x(zci(ii):zci(ii)+1),0)];
            ss = [ss; sign(interp1(y(zci(ii):zci(ii)+1),f(zci(ii):zci(ii)+1),0))];

            % If this is a meshgrid operation, add enough rs
            if numel(r) > 1
                rs = [rs; r(zci(ii))];
            % Otherwise add just one r
            else
                rs = [rs; r];
            end
        end     
    end
end