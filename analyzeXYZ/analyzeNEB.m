#This program is for analysis of trajectory files from Orca (structure optimizations, molecular dynamics, NEB files or relaxed surface scans)
#Programmed by N. Neuman, June 2022
close all

format long

fileName = 'iPrNH2_Ac2O_DMF_qmmm_ZNEBTS_ad_MEP_ALL_trj.xyz';

[Energies,Geometries,nAtoms,kGeom] = readTrajectoryXYZ(fileName);

DistanceVector = zeros(kGeom,1);

#Number of images
nIm = 10;
LastNEB = 290;
nImZOOM = 15;

EnergyMatrix = reshape(Energies(1:LastNEB),nIm,LastNEB/nIm);
EnergyMatrixZOOM = reshape(Energies(LastNEB+1:end),nImZOOM,(kGeom-LastNEB)/nImZOOM);

nNEB = LastNEB/nIm;
nNEBZOOM = (kGeom-LastNEB)/nImZOOM;

figure;
for k = 1:nNEB
    plot((1:nIm)',EnergyMatrix(:,k),'-','Linewidth',1,'Color',[(k/nNEB) (k/nNEB) (k/nNEB)]); hold on
end
for p = 1:nNEBZOOM
    plot([1;nIm/nImZOOM*(2:nImZOOM)'],EnergyMatrixZOOM(:,p),'-','Linewidth',1,'Color',[0.8 (p+1)/(nNEBZOOM+1)/2 (p+1)/(nNEBZOOM+1)/2]); hold on
end
hold off;
ylabel('Energy (Ha)');
xlabel('Image Number');
set(gca, "fontsize", 20);

figure;
surf((1:nIm)'*ones(1,nNEB),(ones(nIm,1)*(nNEB:-1:1)),EnergyMatrix);
xlabel('Image Number');
ylabel('Pass');

figure;
surf((1:nImZOOM)'*ones(1,nNEBZOOM),(ones(nImZOOM,1)*(nNEBZOOM:-1:1)),EnergyMatrixZOOM);
xlabel('Image Number');
ylabel('Pass');
