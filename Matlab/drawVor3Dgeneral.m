function [ I, Avor ] = drawVor3Dgeneral( A1, p, fov, ko, ke )
% Builds the super-res image with help of Voronoi triangulation. The
% density at each data point is estimated as an inverse of the volume of the
% corresponding Voronoi cell. This value is added to the 10-th column of
% the original eventlist yelding 'Avor'. 
% The values of the density are then interpolated to an regular grid with a
% step size of p (nm) using the method 'natural'.
% Optional parameters: fov, ko, ke
if ~exist('ko', 'var')
    ko = 1;
end
if ~exist('ke', 'var')
    ke = size(A1, 1);
end
if ~exist('fov', 'var')
    fov = [max(A1(:,4)) max(A1(:,5)) max(A1(:,6))];
end
if size(A1,2) < 10
    Vor = VorArea3D( A1, ko, ke );
    Area = Vor{1}; %loc.dens
    %A1 = Vor{4};
    ic = Vor{5};
    Areacorr = zeros(size(A1,1),1);
    Areacorr(:) = Area(ic,:);
    Avor = [A1, Areacorr];
elseif size(A1,2) == 10
    Avor = A1;
end
%Avor(:,10) = Avor(:,10) .* Avor(:,7) ./ mean(Avor(:,7)); %multiply by number of photons
% drawing
if exist ('p', 'var')
s = round(fov / p); % size in pixels
M = zeros(prod(s), 3); %template for query matrix
a = 1;
for i = 1:s(1)
    for j = 1:s(2)
        for k = 1:s(3)
            M(a,1) = i;
            M(a,2) = j;
            M(a,3) = k;
            a = a+1;
        end
    end
end
M = M * p; %convert to nanometers

% interpolate Anew to M
M(:,4) = griddata(Avor(:,4), Avor(:,5), Avor(:,6), Avor(:,10), M(:,1), M(:,2), M(:,3), 'natural');

%delete negative values
M(M(:,4) < 0 , 4) = 0;
M(isnan(M(:,4)), 4) = 0;
M(isinf(M(:,4)), 4) = 0;

% build an image from M(x,y,z)
I = zeros(s);
for k = 1 : size(M,1)
    I(M(k,2) / p, M(k,1) / p, M(k,3) / p) = M(k,4);
end
else
    I = [];
end
end
%figure; imshow(10*uint8(256*I/max(max(I)))); colormap(hot);

