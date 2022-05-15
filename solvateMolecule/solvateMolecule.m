clear, close all

#This line generates the cell arrays and matrices used by the solvent generating routine.
#TODO: Later add the number of molecules and size of the solvation box.
[solventMolCell,solute,sizeBox,nSolvMol] = solvation_input();

myMol = solute;

#These variables should be calculated inside solvation_input
xMolMin = min(myMol(:,2));
xMolMax = max(myMol(:,2));
yMolMin = min(myMol(:,3));
yMolMax = max(myMol(:,3));
zMolMin = min(myMol(:,4));
zMolMax = max(myMol(:,4));

xMolCenter = (xMolMax+xMolMin)/2;
yMolCenter = (yMolMax+yMolMin)/2;
zMolCenter = (zMolMax+zMolMin)/2;

#sizeBox = 20;

xBoxMin = -sizeBox/2+xMolCenter;
xBoxMax = +sizeBox/2+xMolCenter;
yBoxMin = -sizeBox/2+yMolCenter;
yBoxMax = +sizeBox/2+yMolCenter;
zBoxMin = -sizeBox/2+zMolCenter;
zBoxMax = +sizeBox/2+zMolCenter;

#nSolvMol = 400;

if iscell(solventMolCell)

  nSolventType = size(solventMolCell,1);

  if nSolventType == 1
    solventMol = solventMolCell{1,1};
    solventMolRatio = solventMolCell{1,2};
  end


end

#Minimal distance allowed between any two atoms, each in different molecules.
#Eventually this will be replaced by a more complex criterion depending on the types of atoms
vanDerWaals = 1.65;

solvatedMolecule = zeros(100000,4); #up to 100000 atoms

#I start adding water molecules at the end of myMol
tNext = size(myMol,1);
solvatedMolecule(1:tNext,:) = myMol;
tIter = 1;
for k = 1:nSolvMol

  currPos = rand(3,1).*[xBoxMax-xBoxMin;yBoxMax-yBoxMin;zBoxMax-zBoxMin]+[xBoxMin;yBoxMin;zBoxMin];
  alpha = rand(1)*pi;
  beta = rand(1)*pi;
  gamma = rand(1)*pi;

  #Function erot can be obtained by installing easyspin.
  RotMat = erot(alpha,beta,gamma);

  if nSolventType == 1
      currentSolventMol = solventMol;
      currentSolventMol = [currentSolventMol(:,1), (RotMat*currentSolventMol(:,2:4)')'];
  elseif nSolventType == 2

      chanceFirst = (solventMolCell{1,2}/(solventMolCell{1,2}+solventMolCell{2,2}));    #Normalize the chance
      dice = rand(1);

      if dice < chanceFirst
          currentSolventMol = solventMolCell{1,1};
      else
          currentSolventMol = solventMolCell{2,1};
      end
      currentSolventMol = [currentSolventMol(:,1), (RotMat*currentSolventMol(:,2:4)')'];
  end

  tStep = size(currentSolventMol,1);

  #Only add a molecule if it doesn't collide with previous ones.
#  if calculateDistance(currPos,currentSolventMol,solvatedMolecule,tNext)>vanDerWaals
  if areTouching(currPos,currentSolventMol,solvatedMolecule,tNext) == 0
    solvatedMolecule(tNext+1:tNext+tStep,:) = currentSolventMol+[0,currPos'];
    tNext = tNext+tStep;
    k = k + 1;
  endif

  tIter = tIter + 1;

  #This is to avoid an infinite loop
  #100 is supposed to be sufficiently large so that all molecules find their place
  if tIter > 100*nSolvMol
    break
  endif

end

solvatedMolecule = solvatedMolecule(1:tNext,:);

#Remember to change the name of your output .xyz file so not to overwrite the previous ones.
#Uncomment next line for saving
#save solvatedMoleculeB.xyz solvatedMolecule -ascii;






