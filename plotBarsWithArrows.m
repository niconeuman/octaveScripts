function [xAxesBars,EnergyBars] = plotBarsWithArrows(Energies,Occupations)
  
  nPoints = 50;
  xvec = linspace(2,3,nPoints);  %A number which is divisible by 2, 3, 4 and 5, to be able to make shorter vectors
  nLevels = length(Energies);
  DeltaEMax = max(Energies)-min(Energies);
  
  deltaThr = DeltaEMax/20; %This is a threshold, below which orbital energies are considered to be too close.
  EnergyBars = zeros(nLevels,nPoints); %At most I will have nLevels, but if there is degeneracy, I will use less
  xAxesBars = zeros(nLevels,nPoints);
  
  arrowScaleY = DeltaEMax/300;
  xArrows = zeros(nLevels,60); %60 is double the points I used for the arrow
  yArrows = zeros(nLevels,60); %60 is double the points I used for the arrow
  currLevel = 1;
  k = 1;
  
  pos = "left";
  
  for t = 1:nLevels-1
    
    l = k + 1;
    #disp('k is now:');
    #disp(k);
        if abs(Energies(l)-Energies(k))<deltaThr
            
            xvec1 = linspace(1.25,2.25,nPoints);
            xvec2 = linspace(2.75,3.75,nPoints);
            
            EnergyBars(currLevel,:) = [Energies(k)*ones(1,nPoints)];
            xAxesBars(currLevel,:) = xvec1;
            
            [xout,yout] = arrow(1.75,Energies(k),Occupations(k),arrowScaleY,pos);
            
            xArrows(currLevel,:) = xout;
            yArrows(currLevel,:) = yout;
            
            currLevel++;
            EnergyBars(currLevel,:) = [Energies(l)*ones(1,nPoints)];
            xAxesBars(currLevel,:) = xvec2;
            
            [xout,yout] = arrow(3.25,Energies(l),Occupations(l),arrowScaleY,pos);
            
            xArrows(currLevel,:) = xout;
            yArrows(currLevel,:) = yout;
            
            currLevel++;
            k = k+1;
            
            if strcmp(pos,"left")
              pos = "right";
            else
              pos = "left";
            endif
        else
            EnergyBars(currLevel,:) = [Energies(k)*ones(1,nPoints)];
            xAxesBars(currLevel,:) = xvec;
            
            [xout,yout] = arrow(2.5,Energies(k),Occupations(k),arrowScaleY,pos);
            
            xArrows(currLevel,:) = xout;
            yArrows(currLevel,:) = yout;
            
            currLevel++;
##            EnergyBars(currLevel,:) = [Energies(l)*ones(1,nPoints)];
##            xAxesBars(currLevel,:) = xvec;
            if strcmp(pos,"left")
              pos = "right";
            else
              pos = "left";
            endif                  
        end
     k = k+1;  %I decouple this index from the looping index so I can increase it from within the loop
     if (k == nLevels)
            EnergyBars(currLevel,:) = [Energies(k)*ones(1,nPoints)];
            xAxesBars(currLevel,:) = xvec;
            
            [xout,yout] = arrow(2.5,Energies(k),Occupations(k),arrowScaleY,pos);
            
            xArrows(currLevel,:) = xout;
            yArrows(currLevel,:) = yout;
            
            if strcmp(pos,"left")
              pos = "right";
            else
              pos = "left";
            endif
       break;
     endif
  endfor
  
  function [xout,yout] = arrow(x,y,Occupation,arrowScaleY)
    
    if strcmp(pos,"left")
      moveX = -0.15;
    else
      moveX = +0.15;
    endif
    
    ScaleX = 0.25;
    ShiftX = 0.03;
    xa = [-0.45 -0.42 -0.39 -0.36 -0.33 -0.30 -0.27 -0.24 -0.21 -0.18   -0.16 -0.16 -0.16 -0.16 -0.16 -0.16 -0.16 -0.16 -0.16 -0.16 -0.16 -0.16 -0.16 -0.16 -0.16 -0.16 -0.16 -0.16 -0.16 -0.16 ];
    ya = [4 5.5 7 8.5 10 11.5 13 14.5 16 17.5  19 17 15 13 11 9 7 5 3 1 -1 -3 -5 -7 -9 -11 -13 -15 -17 -19];
     
    #xa = xa*Scale+x;
    #ya = ya*Scale+y;
    
    if Occupation == 1
      xout = kron(xa,[1 1])*ScaleX+x-ShiftX+moveX; %I double the points so that in each case there are 30 points in total.
      yout = kron(ya,[1 1])*arrowScaleY+y;
    elseif Occupation == -1
      xout = kron(-xa,[1 1])*ScaleX+x-ShiftX+moveX; %I double the points so that in each case there are 30 points in total.
      yout = kron(-ya,[1 1])*arrowScaleY+y;
    elseif Occupation == 2
      xout = [xa*ScaleX-ShiftX -xa*ScaleX+ShiftX]+x+moveX;
      yout = [ya -ya]*arrowScaleY+y;
    else
      xout = 1000*ones(1,60);  #This will be out of scale
      yout = ones(1,60);
    endif
    
  endfunction
  
  
  figure('position',[200,200,600,800]);
  plot(xAxesBars',EnergyBars','sk','MarkerFaceColor','k','MarkerSize',3); hold on;
  #disp(xArrows);
  #disp(yArrows);
  plot(xArrows',yArrows','sr','MarkerFaceColor','r','MarkerSize',2);
  axis([0;5;1.25*min(Energies);1.25*max(Energies)],"tic[y]")
  #axis("off")
  set(gca, "fontsize", 20);
  set(gca, 'XColor', [1 1 1])
  
  box off;
  
  
  
  
end
