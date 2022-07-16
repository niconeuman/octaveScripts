function [Energies,Geometries,nAtoms,kGeom] = readTrajectoryXYZ(fileName)

fileID = fopen(fileName, "r");

#This program can handle up to 2000 structures. This value is truncated to the actual number of structures at the end.
maxGeom = 2000;
Energies = zeros(maxGeom,1);
Geometries = cell(maxGeom,2);

nAtoms = textscan(fileID,"%f",1,"Delimiter","\t");
nAtoms = nAtoms{1};

kGeom = 1;

while (kGeom < maxGeom) && ~feof(fileID)

    nAtomsDiscard = textscan(fileID,"%f",1,"Delimiter","\t");

    Energy = textscan(fileID,"%s %s %s %s %s %f",1);

        if feof(fileID)
            break;
        end

    Energies(kGeom) = Energy{1,6};

    currentGeometry = textscan(fileID,"%2s %f %f %f",nAtoms,"Delimiter","\t");

    currentAtomList = currentGeometry{1};
    currentCoordinates = [currentGeometry{2} currentGeometry{3} currentGeometry{4}];

    Geometries{kGeom,1} = cell(nAtoms,1);

    Geometries{kGeom,1} = currentAtomList;
    Geometries(kGeom,2) = currentCoordinates;

    kGeom += 1;

end

kGeom -= 1;

Energies = Energies(1:kGeom);
Geometries = Geometries(1:kGeom,:);



fclose(fileID);







end

