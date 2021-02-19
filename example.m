
close all

tic
    bifurcation1D(@p02_00);
toc

tic
    bifurcation1D(@p02_01);
toc


tic
    bifurcation1D(@p02_02);
toc

return


function [dxs, fPrimes] = p02_00(xs,r)

    % P1
    
    dxs     = r*xs - xs^2*(1-xs);
    fPrimes = r - 2*xs^1 + 3*xs^2;

end

function [dxs, fPrimes] = p02_01(xs,r)

    % P1
    
    dxs     = r*xs - xs.^2.*(1-xs);
    fPrimes = r - 2.*xs.^1 + 3.*xs.^2;

end


function [dxs, fPrimes] = p02_02(xs,r)

    % P1
    
    dxs     = r.*xs - xs.^2.*(1-xs);
    fPrimes = r - 2.*xs.^1 + 3.*xs.^2;

end