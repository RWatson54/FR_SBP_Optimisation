function [refPoly, nSoln, nOutr, xOutr, symMatx] = getPolygonSymmetry(nSide, symOrb)

%GETPOLYSYMMETRY generates the polygon and sets up the symmetry matrix

% -------------------------------------------------------------

% -- The inputs are:

%       nSide -- The number of sides of the polygon we're working with
%       symOrb -- the count of the different symmetry orbits being used

% -- The outputs are:

%       refPoly -- the polygon structure we're using
%       nSoln -- the number of solution points involved
%       nOutr -- the number of outer points defining the shape involved
%       xOutr -- the coordinates defining the shape
%       symMatx -- the matrix which applies symmetry to the generated points

% -------------------------------------------------------------

% -- Set up the polygon
refPoly = nsidedpoly(nSide,'Center',[0 0],'SideLength',2);

% -- Get the number of solution points
nSoln = symOrb(1) + nSide*symOrb(2) + 2*nSide*symOrb(3);

% -- Number the points anticlockwise rather than clockwise
nOutr = size(refPoly.Vertices,1);
xOutr = flipud(refPoly.Vertices);

% -- Generate the matrix which "does symmetry" automatically over the polygon
RefsymMat = zeros( nOutr*2*2, 2);
for i = 1:nOutr
   iA = i*2-1;
   iB = i*2;
   
   T1 = [cos(2*pi*(i-1)/nOutr), -sin(2*pi*(i-1)/nOutr); sin(2*pi*(i-1)/nOutr),  cos(2*pi*(i-1)/nOutr)];
   
   RefsymMat(iA,:) = T1(1,:);
   RefsymMat(iB,:) = T1(2,:);
   
end
for i = 1:nOutr
   iA = i*2-1 + 2*nOutr;
   iB = i*2   + 2*nOutr;
   
   twot = 2 * atan2(xOutr(1,2),xOutr(1,1));
   
   
   T0 = [cos(twot), sin(twot); sin(twot), -cos(twot)];
   T1 = T0 * [cos(2*pi*(i-1)/nOutr), -sin(2*pi*(i-1)/nOutr); sin(2*pi*(i-1)/nOutr),  cos(2*pi*(i-1)/nOutr)];
   
   RefsymMat(iA,:) = T1(1,:);
   RefsymMat(iB,:) = T1(2,:);
   
end
RefsymMat(:,3) = zeros(size(nOutr*2*2, 1));

% -- Refer the symmetric set to the unit triangle
TranssymMat = [[0 0 1]', [0.5*(xOutr(end,1)+xOutr(1,1)), 0.5*(xOutr(end,2)+xOutr(1,2)), 1]', [xOutr(1,1), xOutr(1,2), 1]'] * inv([[0 0 1]', [1 0 1]', [0 1 1]']);
symMatx = RefsymMat * TranssymMat;

end

