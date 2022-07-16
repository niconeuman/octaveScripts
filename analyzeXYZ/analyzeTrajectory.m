#This program is for analysis of trajectory files from Orca (structure optimizations, molecular dynamics, NEB files or relaxed surface scans)
#Programmed by N. Neuman, June 2022

format long

#fileName = 'CoVID4MC2MCcor_Py_2_H2O_10_MeCN_8_PBE_scan_eh_trj.xyz';
##fileName = 'iPrNH2_Ac2O_DMF_qmmm_scan_ac2_trj.xyz';
##fileName = 'PEG_458_H2O_DCM_B973c_opt_ba_trj.xyz';
fileName = 'iPrNH2_Ac2O_DMF_qmmm_ZNEBTS_ad_MEP_ALL_trj.xyz';

[Energies,Geometries,nAtoms,kGeom] = readTrajectoryXYZ(fileName);

DistanceVector = zeros(kGeom,1);

atom1 = 3;
atom2 = 14;

for k = 1:kGeom
  Pos1 = Geometries{k,2}(atom1,:);
  Pos2 = Geometries{k,2}(atom2,:);

  DistanceVector(k) = norm(Pos2-Pos1);

end

figure;
plot(DistanceVector,Energies,'sk','Linewidth',2);
ylabel('Energy (Ha)');
xlabel('Distance (A)');
set(gca, "fontsize", 20);
