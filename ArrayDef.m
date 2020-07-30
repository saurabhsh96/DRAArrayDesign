%Antenna Array making required TSV
clear;
%% Input definations
%Number of elements in array
Nx = 6;
Ny = 6;
Nz = 0;

%Initial positions of individual elements
xPos = zeros(Nx, 1);
yPos = zeros(Nx, 1);
Z = zeros(Nx*Ny, 1);

%Phase, Amplitude and element vector
%Phase = zeros(Nx*Ny, 1);
Magnitude = ones(Nx*Ny, 1);
Elements = (1:Nx*Ny)';

%Spaceshift calculations
f = 6.5e9;
c = 3e8;
lam = c/f;
k0 = 2*pi/lam;
spaceshift = lam/2;

%Position mapping
startPosX = -spaceshift*(Nx/2-1) - spaceshift/2;
startPosY = -spaceshift*(Ny/2-1) - spaceshift/2;
xPos(1) = startPosX;
yPos(1) = startPosY;

%Nx and Ny are same, Rect array so,
for ind = 2:Nx
    xPos(ind) = xPos(ind-1) + spaceshift;
    yPos(ind) = yPos(ind-1) + spaceshift;
end

%Main matrix
X = zeros(Nx*Ny, 1);
Y = zeros(Nx*Ny, 1);

%Appending elements in Nx*Ny array
k = 1;
for ind = 1:Nx:Nx*Ny
    for indj = 1:Nx
        Y(ind+indj-1) = yPos(k);
        X(ind+indj-1) = xPos(indj);
    end
    k = k+1;
end

%Window calculations
%Gives 25.1 deg beamwidth as required
%winX = chebwin(Nx, 75);
%winY = chebwin(Ny, 75);

%Works for 31 deg beamwidth
winX = chebwin(Nx, 80);
winY = chebwin(Ny, 80);

win = winX*winY';

%Magnitude Tapering
Magnitude = Magnitude.*win(:);

%Angle definitions, I guess the orientation of tha array
Phi = zeros(Nx.*Ny, 1);
Theta = zeros(Nx.*Ny, 1);
Gamma = zeros(Nx.*Ny, 1);

%Template
lines = ["# Created by Saurabh Nerkar\n" ... 
"# On Saturday, April 11, 2020 at 9:14:35 PM\n" ...
"# unit: meters\n" ...
"# design frequency: 6500000000 Hz\n" ...
"# Element\tX\tY\tZ\tMagnitude\tPhase\tPhi\tTheta\tGamma\n"];

%Writing template
fid = fopen('tableData.txt','wt');
for ind = lines
    fprintf(fid, ind);
end

%Phase shift calculations Phi = 0; E-plane; phi = 90 H-plane
%Required beam, phi = 45 deg; theta = 60 deg
drad = pi/180;
theta0 = 60*drad; %page no 352 Balanis
phi0 = 45*drad;
%In this case, spaceshift is same for x and y, defined as above, for cases
%with different spaceshift we can define variables accordingly.
Bx = -k0.*spaceshift.*sin(theta0).*cos(phi0)./drad;
By = -k0.*spaceshift.*sin(theta0).*sin(phi0)./drad;

Phase = zeros(Nx, Ny);
for indy = 1:Nx
    for ind = 1:Ny
        Phase(ind, indy) = ind.*Bx + indy.*By;
        reqInt = floor(Phase(ind, indy)/360);
        Phase(ind, indy) = Phase(ind, indy) - 360*reqInt;
    end
end

T = [Elements X Y Z Magnitude Phase(:) Phi Theta Gamma];
%Writing data matrix
fmt = '%f\t%2.8f\t%2.8f\t%2.8f\t%2.8f\t%2.8f\t%2.8f\t%2.8f\t%2.8f\n';
for ind = 1:size(T, 1)
    %fprintf(fid, [num2str(T(ind,:)),'\n']);
    fprintf(fid, fmt, T(ind,:));
end
fclose(fid);