function centerOfMass = centerOfMass(solventMol);

  #Take a molecule with array structure:

  #mol = [AN1 x1 y1 z1
  #       AN2 x2 y2 z2
  #       ...
  #       ...
  #       ANN xN yN zN];

  #And calculate its center of Mass

  sumMass = sum(solventMol(:,1));

  centerOfMass = [0 sum(solventMol(:,2).*solventMol(:,1))/sumMass, sum(solventMol(:,3).*solventMol(:,1))/sumMass, sum(solventMol(:,4).*solventMol(:,1))/sumMass];



end


