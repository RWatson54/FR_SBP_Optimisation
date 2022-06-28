function [xFInt, wFInt, xFNrm] = getFaceIntegration(xOutr, nFInt)

% This function computes the integration weights (and face normals) around a polygon defined by xOutr

% -------------------------------------------------------------

% -- The inputs are:

%       xOutr -- the coordinates defining the outer points and shape
%       nFInt -- the number of points to use on each face

% -- The outputs are:

%       xFInt -- the coordinates of the integration points along each face
%       wFInt -- the weights of the integration points along each face
%       xFNrm -- the face normals for each of the integration points

% -------------------------------------------------------------

% -- Get the number of faces from the outer points
nFace = size(xOutr,1);
for iFace = 1:nFace

    % -- Work out which of the outer points bound this face
    lP = mod(iFace  -1,nFace)+1;
    rP = mod(iFace+1-1,nFace)+1;

    % -- Get the direction of the normal vector
    n(1) = (xOutr(rP,2) - xOutr(lP,2));
    n(2) = (xOutr(lP,1) - xOutr(rP,1));

    % -- Get the integration weights along the edge
    [xFInt((iFace-1)*nFInt+1:iFace*nFInt,:), wFInt((iFace-1)*nFInt+1:iFace*nFInt,:)] = lgwt(nFInt, xOutr(lP,:), xOutr(rP,:));

    % -- And add in the normalised normals
    xFNrm((iFace-1)*nFInt+1:iFace*nFInt,1) = n(1) / norm(n);
    xFNrm((iFace-1)*nFInt+1:iFace*nFInt,2) = n(2) / norm(n);

end

end

