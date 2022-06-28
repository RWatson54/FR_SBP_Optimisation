% -- Tidy up
close all
clear all

% -- Set up some counters
Side = [2 3 4 5 6 7];
Maximal2D = [4 9 16 25 26 49];
Total2D = [3 6 10 15 21 28];
Euclid2D = [3 6 11 17 26 35];

% -- Triangles
% ============
nSide = 3;

% -- Loop up to a large number of points, working out which ones have appropriate
possible_nSoln = (1:100);

% -- Get the numbers of points which are capable of giving symmetric solutions
plausible_nSoln = find(mod(possible_nSoln , nSide)<2);

% -- Get the point sets which are compatible for each type of points
[nSoln_Maximal, ~, iB] = intersect(plausible_nSoln, Maximal2D);
nFace_Maximal = Side(iB);
[nSoln_Total, ~, iB] = intersect(plausible_nSoln, Total2D);
nFace_Total = Side(iB);
[nSoln_Euclid, ~, iB] = intersect(plausible_nSoln, Euclid2D);
nFace_Euclid = Side(iB);

% -- Get the available run inputs for each set
nO = 0; nTrial = 144;
% -- -- Maximal Bases
for iS = 1:size(nSoln_Maximal,2)
    Orbs = getPlausibleOrbits(nSoln_Maximal(iS), nSide);
    for iO = 1:size(Orbs,1)
        nO = nO + 1;
        RunInput(nO).nSide = nSide;
        RunInput(nO).nSoln = nSoln_Maximal(iS);
        RunInput(nO).SymOrb = Orbs(iO,:);
        RunInput(nO).nFace = nFace_Maximal(iS);
        RunInput(nO).nTrial = 144;
        RunInput(nO).compBasisType.Type = 'Maximal2D';
        RunInput(nO).compBasisType.Eps  = 'Dummy';
        RunInput(nO).compBasisType.xC   = 'Dummy';
        RunInput(nO).testBasisType.Type = 'Maximal2D';
        RunInput(nO).testBasisType.Eps  = 'Dummy';
        RunInput(nO).testBasisType.xC   = 'Dummy';
        RunInput(nO).wT = ones(nTrial,1);
    end
end
% -- -- Total Bases
for iS = 1:size(nSoln_Total,2)
    Orbs = getPlausibleOrbits(nSoln_Total(iS), nSide);
    for iO = 1:size(Orbs,1)
        nO = nO + 1;
        RunInput(nO).nSide = nSide;
        RunInput(nO).nSoln = nSoln_Total(iS);
        RunInput(nO).SymOrb = Orbs(iO,:);
        RunInput(nO).nFace = nFace_Total(iS);
        RunInput(nO).nTrial = 144;
        RunInput(nO).compBasisType.Type = 'Total2D';
        RunInput(nO).compBasisType.Eps  = 'Dummy';
        RunInput(nO).compBasisType.xC   = 'Dummy';
        RunInput(nO).testBasisType.Type = 'Maximal2D';
        RunInput(nO).testBasisType.Eps  = 'Dummy';
        RunInput(nO).testBasisType.xC   = 'Dummy';
        RunInput(nO).wT = ones(nTrial,1);
    end
end
% -- -- Euclid Bases
for iS = 1:size(nSoln_Euclid,2)
    Orbs = getPlausibleOrbits(nSoln_Euclid(iS), nSide);
    for iO = 1:size(Orbs,1)
        nO = nO + 1;
        RunInput(nO).nSide = nSide;
        RunInput(nO).nSoln = nSoln_Euclid(iS);
        RunInput(nO).SymOrb = Orbs(iO,:);
        RunInput(nO).nFace = nFace_Euclid(iS);
        RunInput(nO).nTrial = 144;
        RunInput(nO).compBasisType.Type = 'Euclid2D';
        RunInput(nO).compBasisType.Eps  = 'Dummy';
        RunInput(nO).compBasisType.xC   = 'Dummy';
        RunInput(nO).testBasisType.Type = 'Maximal2D';
        RunInput(nO).testBasisType.Eps  = 'Dummy';
        RunInput(nO).testBasisType.xC   = 'Dummy';
        RunInput(nO).wT = ones(nTrial,1);
    end
end

% -- Now, loop over each problem in turn
for iR = 1:size(RunInput,2)
    [RunOutput(iR).xT, RunOutput(iR).E] = generateErrorOptimisedPoints(RunInput(iR).nSide, RunInput(iR).SymOrb, RunInput(iR).nFace, RunInput(iR).nTrial, RunInput(iR).compBasisType, RunInput(iR).testBasisType, RunInput(iR).wT, 2, 10000, 500);
end


% % -- Turn each of these into the required symmetry orbits
% RunType(1).nSoln;
% RunType(1).SymOrb;
% RunType(1).nTrial;
% RunType(1).nFace;
% RunType(1).compBasisType.Type = 'Maximal2D';
% RunType(1).compBasisType.Eps  = 'Dummy';
% RunType(1).compBasisType.xC   = 'Dummy';
% RunType(1).testBasisType.Type = 'Maximal2D';
% RunType(1).testBasisType.Eps  = 'Dummy';
% RunType(1).testBasisType.xC   = 'Dummy';
% RunType(1).wT = ones(nTrial,1);

% -- And run the code with these settings