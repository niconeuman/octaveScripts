function are_touching = areTouching(currPos,solventMol,solvatedMolecule,tNext)

  #currPos = [xCofM; yCofM; zCofM];

  Molecule1 = solventMol+[0,currPos'];

  xSizeMol = max(solventMol(:,2))-min(solventMol(:,2));
  ySizeMol = max(solventMol(:,3))-min(solventMol(:,3));
  zSizeMol = max(solventMol(:,4))-min(solventMol(:,4));

  rMol = sqrt(xSizeMol^2+ySizeMol^2+zSizeMol^2);                                #Assuming a spheroidal molecule

  #This is presumably bigger, or much bigger, than the new solvent molecule.
  Molecule2 = solvatedMolecule(1:tNext,:);

  #Define the limits of the solvated region
  xSolvatedMin = min(Molecule2(:,2));
  xSolvatedMax = max(Molecule2(:,2));
  ySolvatedMin = min(Molecule2(:,3));
  ySolvatedMax = max(Molecule2(:,3));
  zSolvatedMin = min(Molecule2(:,4));
  zSolvatedMax = max(Molecule2(:,4));

  length1 = size(Molecule1,1);
  length2 = size(Molecule2,1);

  are_touching = 0; #Very far away!
  maybe_touching = 0; #Very far away!

                      #H  #He  #Li  #Be  #B   #C   #N   #O   #F  #Ne  #Na  #Mg  #Al  #Si  #P   #S   #Cl  #Ar
  vanderWaalsRadii = [1.2 1.0 2.08 2.08 2.08 1.85 1.75 1.70 1.65 1.56 2.25 2.40 2.39 2.32 2.23 2.14 2.06 1.97];

  #If the CofM of the new molecule is inside the region where the solvated molecules are (there could be empty spaces), they may be touching
  if (currPos(1)>xSolvatedMin) && (currPos(1)<xSolvatedMax)
    if (currPos(2)>ySolvatedMin) && (currPos(2)<ySolvatedMax)
      if (currPos(3)>zSolvatedMin) && (currPos(3)<zSolvatedMax)
        maybe_touching = 1;
      endif
    endif
  endif


  if maybe_touching == 1  #If there is an overlap in the region, then check better

    for l = 1:length2 #Run over longer molecule first
      if (Molecule2(l,1) < 19)
          VdW2 = vanderWaalsRadii(Molecule2(l,1));
      elseif (Molecule2(l,1) < 30)
          VdW2 = 2.6;
      else
          VdW2 = 2.9;
      end

      #Calculate the distance of all atoms in solvated system to currPos
      currDistToCofM = norm(Molecule2(l,2:4)-currPos');
        if currDistToCofM < rMol

          for k = 1:length1

            if (Molecule1(k,1) < 19)
                VdW1 = vanderWaalsRadii(Molecule1(k,1));
            elseif (Molecule1(k,1) < 30)
                VdW1 = 2.6;
            else
                VdW1 = 2.9;
            end

            #This runs over all atoms in the new molecule and all atoms in the solvated system. It is very slow
            currDistance = norm(Molecule1(k,2:4)-Molecule2(l,2:4));

            if currDistance < (VdW1+VdW2)/2;
              are_touching = 1;
            endif

          endfor
        else
          are_touching = 0;
        endif

    endfor

  endif




end
