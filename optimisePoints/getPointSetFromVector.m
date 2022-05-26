function [xOut] = getPointSetFromVector(symMat, sOrb, sVec)

%GETPOINTS generates the coordinates of the solution points of a polygon from an input vector

% This might not be a bad target for some optimisation/vectorisation, since it's called frequently

% -------------------------------------------------------------

% -- The inputs are:

%       symMat -- The previous-calculated symmetry matrix for performing the rotations and reflections
%       sOrb -- the vector which enumerates the symmetry orbits
%       sVec -- the vector which defines the positions of the points within the orbits (which is what is ultimately optimised)

% -- The outputs are:

%       xOut -- the coordinates of the symmetrical points on the polygon

% -- More details:

%      sVec is defined across the unit triangle [0;0  ,  0;1   ,    1;0]
%      and returns, for any sVec values, a symmetric set of points

%      There are three symmetry orbits on any polygon:
%          1 - the central point - 1 point (0,0)
%          2 - the points along the lines (0-1,0) and (0,0-1)
%          3 - the general points (0-1,0-1)

% Get the number of sides
nSide = size(symMat,1) / 4;

% Set the sizes of each orbit
OrbLength = [1, nSide, 2*nSide];

% sOrb = say [ 1, 2, 2 ]
xOut = zeros(dot(sOrb,OrbLength),2);

iX = 0; iF = 0;
% Muck around with first symmetry orbit
for iOrb = 1:sOrb(1)
    % Get the memory locations
    iX = iX + 0;
    iS = iF + 1;
    iF = iF + OrbLength(1);
    
    % Straight up just set this to zero
    xOut(iS:iF,:) = [0, 0];
end
% Muck around with second symmetry orbit
for iOrb = 1:sOrb(2)
    % Get the memory locations
    iX = iX + 1;
    iS = iF + 1;
    iF = iF + OrbLength(2);
    
    % Limit the values from -1 to 1 of sVec(iX,1)
    if sVec(iX) < 0
        xL = [min(max(-sVec(iX),0.01),1), 0];
    else
        xL = [0, min(max( sVec(iX),0.01),1)];
    end
    xOut(iS:iF,:) = reshape(symMat(1:nSide*2,  :) * [xL'; ones(1,size(xL,1))],2,[])';
end
% Muck around with third symmetry orbit
for iOrb = 1:sOrb(3)
    % Get the memory locations
    iX = iX + 2;
    iS = iF + 1;
    iF = iF + OrbLength(3);
    
    % Limit the values from sVec to the properly within triangle (0.01, 0.01; 0.99, 0.01; 0.01, 0.99)
    xL = [max(min(sVec(iX-1),0.99),0.01), min(max(sVec(iX),0.01),1.0-max(min(sVec(iX-1),0.99),0.01))];
    xOut(iS:iF,:) = reshape(symMat(1:nSide*2*2,  :) * [min(xL,1)'; ones(1,size(xL,1))],2,[])';
end

end

