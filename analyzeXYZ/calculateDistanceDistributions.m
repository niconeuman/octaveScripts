function [distanceDistributions,typeA,typeB] = calculateDistanceDistributions(Geometries,typeA,typeB)

#This function reads the Geometries cell array and for each structure generates a distribution of the A-B distances

distanceDistributions = cell(size(Geometries,1),1);

if length(typeA) == 1
typeA = [typeA, ' '];
end

if length(typeB) == 1
typeB = [typeB, ' '];
end

##disp(typeA);
##disp(typeB);

atomList = Geometries{1,1};
nAtoms = length(atomList);
maxGeom = size(Geometries,1);

##disp(nAtoms);

#This works. It gives the indices of all cell elements with the matching string
idxListA = find(strcmp(atomList,{typeA}));
idxListB = find(strcmp(atomList,{typeB}));

#If the atom types are different (i.e. 'C' and 'H'), I can loop through idxListA and idxListB without problem.
#But if they are the same atom, I need to loop uniquely (I can't calculate the distance of an atom from itself)

nAtomsA = length(idxListA);
nAtomsB = length(idxListB);

##disp(idxListA);
##disp(idxListB);

pairIndex = [kron(idxListA,ones(size(idxListB))), kron(ones(size(idxListA)),idxListB)];

if (nAtomsA == nAtomsB && all(idxListA == idxListB))
  pairIndex = pairIndex(pairIndex(:,1)~=pairIndex(:,2),:);
end

##disp(pairIndex);

#In principle now I have all unique pairs of indices. I can loop over structures.

distanceList = zeros(size(pairIndex,1),1);

for kGeom = 1

    for kAtom = 1:length(distanceList)
        distance = norm(Geometries{kGeom,2}(pairIndex(kAtom,1),:)-Geometries{kGeom,2}(pairIndex(kAtom,2),:));
        distanceList(kAtom) = distance;


    endfor

    distanceDistributions(kGeom) = distanceList;

end





##for k = 1:nAtomsA-1
##    for l = k+1:nAtomsA
##
##        idxA = idxListA(k);
##        idxB = idx
##
##
##
##    endfor
##
##endfor











end

