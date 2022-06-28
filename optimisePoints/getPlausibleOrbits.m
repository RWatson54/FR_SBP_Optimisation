function [bOrb, bErr] = getPlausibleOrbits(nSoln, nSide)
%GETORBITS 
%   This subroutine gets all the plausible symmetry orbits for the element,
%   based on the number of requested solution points. If it doesn't work,
%   it returns the error value as 1

% -- First, test if this is possible
iR = mod(nSoln,nSide);
if iR > 1
    bErr = 1;
    bOrd = 0;
    return
else
    bErr = 0;
end

% -- Then work out how many possibilities there are
n3 = floor((nSoln-iR) / (2*nSide));

% -- Loop over those possibilities
for iP = 0:n3
    bOrb(iP+1,:) = [iR, (nSoln-2*iP*nSide-iR)/nSide , iP];
end


end

