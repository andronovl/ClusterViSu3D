%% load data in the ascii format, in the 4,5,6 columns are the X,Y,Z coordinates of the fluorophores (in nm);
clear all;
%
A = importdata('H2B.ascii'); %put here your file's name
%% Check the axial distribution of points
[counts, centers] = hist(A(:,6),100);
[muhat,~] = normfit(A(:,6)); %in order to get sigma from the experimental data put here [muhat, sigm] and comment next line
sigm = 300; % standard deviation of the normal distribution which models the dinimition of the number of events while going out of focus
norm = normpdf(centers,muhat,sigm)/(normcdf(581.5,muhat,sigm) - normcdf(-615,muhat,sigm)); %normalize gaussian pdf between -615 and 581.5 - check the Z-extent of your data
plot(centers, counts/(size(A,1)*(centers(2) - centers(1))), centers, norm, 'linewidth', 2);
xlabel('Z, nm','fontsize', 16);
ylabel('pdf','fontsize', 16);
legend('Z proj., \Deltaz = 12 nm', 'ND, \sigma = 300 nm','fontsize', 16); %\Deltaz = centers(2) - centers(1)
title({'Approximation of the Z-axis projection with the normal distribution'});
set(gca,'fontsize', 16, 'LineWidth',1);
%% print if needed, check filename
print('Zprofile.tif', '-dtiff', '-r300');
%% show scattered plot
Data = A(:,4:6);
%
figure; scatter3(Data(:,1), Data(:,2), Data(:,3),'.'); axis equal
%
%print('scatterplot.tif', '-dtiff', '-r300'); %can print it
%%   Find volumes of Voronoi cells
tic;
Data1 = A;
Voronoi = VorArea3D(Data1);
S = 1./Voronoi{1}; %volumes
toc;
%% XY projection image for ROI determination
I = drpar(A, 20); %(dataset, pixelsize)
I = uint8(3*255*I/max(max(I))); %normalize the image as you wish
% Select ROI
g = figure; imshow(I); colormap hot;
h = imfreehand(gca); %can use another selection tool here such as imellipse, etc 
wait(h);
BWm = createMask(h);
close(g);
%% Save the XY projection image if needed, check the filename
imwrite(I, hot(256), 'Imhist.tif');
%% Monte-Carlo simulations - can be slow! a progress bar will appear
tic;
[Histograms, Aroi] = VoronoiMonteCarlo3D( A, BWm, 50, 95, sigm ); %(dataset, mask, number of iterations, confidence level in %, standard deviation of the axial distribution modelled as a ND)
toc;
%%   Plot curves; find the intersections manually and put them as thresh and thresh2 in lines 59,60. 
figure; plot(Histograms(:,1), Histograms(:,2)/(size(Aroi,1)*(Histograms(2,1) - Histograms(1,1))), Histograms(:,1), Histograms(:,4)/(size(Aroi,1)*(Histograms(2,1) - Histograms(1,1))), '--m', Histograms(:,1), Histograms(:,3)/(size(Aroi,1)*(Histograms(2,1) - Histograms(1,1))), '-.g', Histograms(:,1), Histograms(:,5)/(size(Aroi,1)*(Histograms(2,1) - Histograms(1,1))), '--m', 'LineWidth', 2);
set(gca, 'TickLength', [0.02,0.05], 'FontSize', 14);
legend('Data', 'Confidence envelope', 'Mean of randomized distribution');
xlabel('volume, nm^3','FontSize', 16);
ylabel('pdf','FontSize', 16);
title({'Distibution of volumes of Voronoi polyhedrons'},'FontSize', 16);
set(gca,'FontSize', 16, 'LineWidth',1);
%% print the curves if needed, check filename
print('VorCurves.tif', '-dtiff', '-r300');
%   find intersection
%[xi, ~] = polyxpoly(centers, counts, centers, MeanCounts);
%%   find intersections
thresh = 1.085E5; %put here the first intersection between the blue and the green lines
thresh2 = 6.9E5; %put here the second intersection

% comparison with theoretical distributions of Voronoi volumes
% binsize = Histograms(2,1) - Histograms(1,1);
% Vav = (bwarea(BWm) * 20^2 * 1000)/length(Aroi);
% figure; plot(Histograms(:,1), Histograms(:,2) * Vav / (length(Aroi) * binsize), Histograms(:,1), Histograms(:,4) * Vav / (length(Aroi) * binsize), '--m', Histograms(:,1), Histograms(:,3) * Vav / (length(Aroi) * binsize), '-.g', Histograms(:,1), Histograms(:,5) * Vav / (length(Aroi) * binsize), '--m', Histograms(:,1), 3125/24*(Histograms(:,1)/Vav).^4.*exp(-5*Histograms(:,1)/Vav), ':k', 'LineWidth', 2);
% legend('Data', 'Confidence envelope', 'Mean of randomized distribution', '', 'Theoretical curve (lz = 1000 nm)');
% xlabel('volume, nm^3','FontSize', 16);
% ylabel('pdf','FontSize', 16);
% title({'H2B-Alexa647; Leica SR GSD 3D (Wetzlar); 20140226-27','Distibution of volumes of Voronoi polyhedra'},'FontSize', 16);
% set(gca,'FontSize', 16, 'LineWidth',1);
% print('\\DISKSTATION\GSD-Data\Demo GSD 3D Wetzlar 20140226-27\H2B647\5\VorCurvesTheorZ1000.tif', '-dtiff', '-r300');
% %% keep points based on polyhedra
% Acl = Voronoi{4}(S<thresh,:);
% Acl = [zeros(size(Acl,1),3) Acl zeros(size(Acl,1),3)];
% %%
% dlmwrite('VlessThanFixedThresh.ascii', Acl, '\t');
% %%
% Idenoised = drpar(Acl, 20);
% Idenoised = uint8(3*255*Idenoised/max(max(Idenoised)));
% figure; imshow(Idenoised); colormap(hot);%% keep points based on polyhedra
% %%
% imwrite(Idenoised, hot(256), '\\DISKSTATION\GSD-Data\Demo GSD 3D Wetzlar 20140226-27\TubAl647\I_VlessThanFixedThresh.tif');
%% Correct for Z distribution density
mid = norminv((normcdf(0,0,300)+normcdf(600,0,300))/2, 0, 300); %number of events with Z from -mid to mid equals to the number of events with Z from -600 to -mid plus the number of events from mid to 600.
Data = Voronoi{4};
VvorCorr = S.* normpdf(Data(:,3), muhat, 300)/normpdf(mid, 0, 300);
Acl = Data(VvorCorr<thresh,:);
Acl = [zeros(size(Acl,1),3) Acl zeros(size(Acl,1),3)];
%% can save the dataset without "noise" - check filename
dlmwrite('VlessThanCorrectedThresh.ascii', Acl, '\t');
%% image of the denoised dataset
Idenoised = drpar(Acl, 20);
Idenoised = uint8(3*255*Idenoised/max(max(Idenoised)));
figure; imshow(Idenoised); colormap(hot);
%% can save the image, check filename
imwrite(Idenoised, hot(256), 'I_VlessThanCorrectedThresh.tif');
%% show points as a scattered plot in three colors
Colcl = repmat([0 0.5 0], size(Acl,1), 1);
Ard = Voronoi{4}(VvorCorr<thresh2&VvorCorr>thresh,:);
Colrd = repmat([0.85 0.85 0], size(Ard,1), 1);
Abg = Voronoi{4}(VvorCorr>thresh2,:);
Colbg = repmat([0.85 0.85 0.85], size(Abg,1), 1);
Ashow = [Acl(:,4:6); Ard; Abg];
Colshow = [Colcl; Colrd; Colbg];
figure; scatter3(Ashow(:,1), Ashow(:,2), Ashow(:,3),[], Colshow, '.'); axis equal;
%print('filename.tif', '-dtiff', '-r300');
%% Interpolate to an image - can be slow
DensVorCorr = 1./VvorCorr;
Data2 = [zeros(size(Data,1),3) Data zeros(size(Data,1),3) DensVorCorr];
Data2(:,6) = Data2(:,6)+600; % Z coordinates should be positive here (the negative coordinates will be lost)
tic
[IvorInterpGaussCorr, Avor] = drawVor3Dgeneral( Data2, 20, [18000,18000,1200]); %(data, pixelsize, ROI)
toc
%% can save the density map as a .mrc volume (for EM data processing software)
addpath('EMIODist'); % find EMIODist here https://github.com/nikomin/EMIODist
WriteMRC(single(IvorInterpGaussCorr), 200, 'IvorInterpGaussCorr.mrc'); %check file name; 200 = voxel size in Angstroms
%% 3D histogram volume
Ihist = drawHist3D(Data2, 20, [18000,18000,1200]); %(data, voxelsize, ROIxyz)
%% Binarize the 3D density map
InvThresh = 1/thresh;
BW = IvorInterpGaussCorr;
BW = BW >= InvThresh;
figure; patch(isosurface(IvorInterpGaussCorr,InvThresh), 'FaceColor', [0.2 0.5 0.2],'FaceAlpha',.8,'LineStyle', 'none'); axis equal
%% save the image of the segmented volume
print('filename1234.tif', '-dtiff', '-r300');
%% find properties of clusters
stats = regionprops(BW, 'Area', 'Centroid', 'BoundingBox', 'SubarrayIdx', 'Image', 'FilledImage', 'FilledArea', 'PixelIdxList', 'PixelList');
Volumes = cat(1, stats.Area);
Volumes = 8000*Volumes; %8000 = the voxel volume in nm^3
mean(((3*Volumes)/(4*pi)).^(1/3)) %mean diameter of the clusters
std(((3*Volumes)/(4*pi)).^(1/3)) %standard deviation of the radius of the clusters
%patch(isosurface(I,InvThresh))
%build a 3D histogram (size of BW) and count events from PixelIdxList
%% Number of localizations per cluster
for i = 1:size(stats,1)
PixelIdxList = stats(i).PixelIdxList;
N(i) = sum(Ihist(PixelIdxList));
end
mean(N)
std(N)
%% Save for ViSP, intensity = loc.density*10^9
Data3 = [zeros(size(Data,1),3) Data zeros(size(Data,1),3) DensVorCorr];
B = zeros(size(Data3,1), 5);
B(:,1:3) = Data3(:,4:6);
B(:,4) = Data3(:,10)*1000e6;
dlmwrite('LocDensGaussCorr.3d', B, '\t');
%%




%% diameters
hist(2*((3*Volumes)/(4*pi)).^(1/3), 150); % diameters
xlim([25 250]);
ylim([25 3600]);
xlabel('diameter, nm','FontSize', 16);
ylabel('frequency','FontSize', 16);
title({'diameter of H2B clusters'},'FontSize', 14);
set(gca,'FontSize', 16, 'LineWidth',1);
%print('HistVolumes.tif', '-dtiff', '-r300');
%
median(2*((3*Volumes)/(4*pi)).^(1/3))
%% nb of loc-s
hist(N(N<50),50);
xlim([0 50]);
ylim([0 2400]);
xlabel('number of localizations, nm','FontSize', 16);
ylabel('frequency','FontSize', 16);
title({'number of localizations per H2B cluster'},'FontSize', 14);
set(gca,'FontSize', 16, 'LineWidth',1);
%print('HistNbOfEvents.tif', '-dtiff', '-r300');