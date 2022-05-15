function are_touching = areTouching(currPos,solventMol,solvatedMolecule,tNext)
  
  
  Molecule1 = solventMol+[0,currPos'];
  
  Molecule2 = solvatedMolecule(1:tNext,:);
  
  length1 = size(Molecule1,1);
  length2 = size(Molecule2,1);
  
  are_touching = 0; #Very far away!
  
                      #H  #He  #Li  #Be  #B   #C   #N   #O    #F
  vanderWaalsRadii = [1.2 1.0 2.08 2.08 2.08 1.85 1.75 1.70 1.65];
  
  for k = 1:length1
    
    if (Molecule1(k,1) < 10)
        VdW1 = vanderWaalsRadii(Molecule1(k,1));
    elseif (Molecule1(k,1) < 20)
        VdW1 = 2.05;
    elseif (Molecule1(k,1) < 30)
        VdW1 = 2.6;
    else
        VdW1 = 2.9;
    end 
    for l = 1:length2
      if (Molecule2(l,1) < 10)
          VdW2 = vanderWaalsRadii(Molecule1(k,1));
      elseif (Molecule2(l,1) < 20)
          VdW2 = 2.05;
      elseif (Molecule2(l,1) < 30)
          VdW2 = 2.6;
      else
          VdW2 = 2.9;
      end       
      
##      if ((Molecule1(k,1) == 1) && (Molecule2(l,1) == 1))
##        vanDerWaals = 1.2;
##      elseif ((Molecule1(k,1) > 1) && (Molecule2(l,1) == 1)) || ((Molecule1(k,1) == 1) && (Molecule2(l,1) > 1))
##        vanDerWaals = 1.8;
##      elseif (Molecule1(k,1) > 1) && (Molecule2(l,1) > 1)
##        vanDerWaals = 2.4;
##      elseif ((Molecule1(k,1) > 10) && (Molecule2(l,1) <= 10)) || ((Molecule1(k,1) <= 10) && (Molecule2(l,1) > 10))
##        vanDerWaals = 2.8;
##      elseif (Molecule1(k,1) > 10) && (Molecule2(l,1) > 10)
##        vanDerWaals = 3.2;
##      endif
      
      currDistance = norm(Molecule1(k,2:4)-Molecule2(l,2:4));
      
      if currDistance < (VdW1+VdW2)/2;
        
        are_touching = 1;
        
      endif
      
    endfor
  endfor
  
  
  
  
  
end