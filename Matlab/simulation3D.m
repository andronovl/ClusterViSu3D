%% set values FOV = 1000*1000*1000 nm
clear all;
r = 30; % radius of clusters (nm)
bkgdDens = 0.6E-4; % background density (nm^-3)
clstDens = 4E-4; % density in clusters (nm^-3)
Nbclusters = 175; % Number of clusters
%% generate (original) data - regularly distributed clusters
% background = rand(round(bkgdDens * (1000)^3),3) * (1E3); 
% [centersX, centersY, centersZ] = meshgrid(100:200:900, 100:200:900, 100:200:900);
% centers = reshape(cat(3, centersX, centersY, centersZ), 125, 3)-50;
% cluster=[];
% for i = 1:size(centers,1) % for each cluster
% clusterPol = rand(round(8 * r^3 * (clstDens - bkgdDens)),3) * r * 2 ;
% clusterPol1 = clusterPol;
% c = 0;
% for j = 1:size(clusterPol1,1)
%     if (clusterPol1(j,1)-r)^2 + (clusterPol1(j,2)-r)^2 + (clusterPol1(j,3)-r)^2 > r^2
%         clusterPol(j-c,:) = [];
%         c = c + 1;
%     end
% end
% clusterXY = [centers(i,1) + clusterPol(:,1), centers(i,2) + clusterPol(:,2), centers(i,3) + clusterPol(:,3)];
% cluster = [cluster; clusterXY];
% end
% Data = [background; cluster];
%% generate (original) data - randomly distributed clusters with a constrain: distance between the centers of clusters >= 5r; distance from the borders to the centers of the clusters >= 2.5r
% background = rand(round(bkgdDens * (1000)^3),3) * (1E3); 
% %centers = rand(125,3)*800+50;
% centers = zeros(175,3);
% cluster=[];
% for i = 1:size(centers,1) % for each cluster
%     centers(i,:) = rand(1,3)*(1000-5*r)+1.5*r;
%     while min(min(pdist(centers(1:i,:)))) < 5 * r
%         centers(i,:) = rand(1,3)*(1000-5*r)+1.5*r;
%     end        
% clusterPol = rand(round(8 * r^3 * (clstDens - bkgdDens)),3) * r * 2 ;
% clusterPol1 = clusterPol;
% c = 0;
% for j = 1:size(clusterPol1,1)
%     R = sqrt((clusterPol1(j,1)-r)^2 + (clusterPol1(j,2)-r)^2 + (clusterPol1(j,3)-r)^2);
%     if R > r
%         clusterPol(j-c,:) = [];
%         c = c + 1;
%     end
% end
% clusterXY = [centers(i,1) + clusterPol(:,1), centers(i,2) + clusterPol(:,2), centers(i,3) + clusterPol(:,3)];
% cluster = [cluster; clusterXY];
% end
% Data = [background; cluster];
%% generate (original) data - randomly distributed clusters with a constrain: distance between the centers of clusters >= 5r; distance from the borders to the centers of the clusters >= 2.5r; smooth borders (from erf function); generate exact number of points in clusters (-/+ background)
background = rand(round(bkgdDens * (1000)^3),3) * (1E3); 
%centers = rand(125,3)*800+50;
centers = zeros(Nbclusters,3);
cluster=[];
for i = 1:size(centers,1) % for each cluster
    centers(i,:) = rand(1,3)*(1000-5*r)+2.5*r;
    while min(min(pdist(centers(1:i,:)))) < 5 * r
        centers(i,:) = rand(1,3)*(1000-5*r)+2.5*r;
    end
PtsInClust = round(4/3*pi*r^3*(clstDens-bkgdDens) + rand - 0.5);
elevation = asin(2*rand(PtsInClust,1)-1);
azimuth = 2*pi*rand(PtsInClust,1);
radii = [];
while size(radii,1) < PtsInClust
    R = 2*r*(rand.^(1/3));
    probtokeep = (erf(4*(r-R)/r)+1)/2;
    if rand < probtokeep
        radii = [radii; R];
    end
end
clusterPol = zeros(PtsInClust,3);
[clusterPol(:,1),clusterPol(:,2),clusterPol(:,3)] = sph2cart(azimuth,elevation,radii);
clusterXY = [centers(i,1) + clusterPol(:,1), centers(i,2) + clusterPol(:,2), centers(i,3) + clusterPol(:,3)];
cluster = [cluster; clusterXY];
end
Data = [background; cluster];
%% scatter plot of the data
figure; scatter3(Data(:,1), Data(:,2), Data(:,3),15,'.'); axis equal
view(20,10);
set(gca, 'TickLength', [0.01,0.01], 'FontSize', 14, 'LineWidth', 0.5, 'Xtick', 0:200:1000, 'Ytick', 0:500:1000, 'Ztick', 200:200:1000);
xlabel('nm','FontSize', 14);
ylabel('nm','FontSize', 14);
zlabel('nm','FontSize', 14);
%title({'R_c_l=30nm; N=175(cl); d_c_l=3\times10^-^4; d_{bgd}=0.6\times10^-^4','Distibution of volumes of Voronoi polyhedrons'},'FontSize', 16);
xlim([0 1000]);
ylim([0 1000]);
zlim([0 1000]);
%% print figure if needed, check the name
print('r30bg0.6E-4cl4E-4N175erf.tif', '-dtiff', '-r400');
%%   Find volumes of Voronoi cells, Monte-Carlo
Data1 = zeros(size(Data,1),9);
Data1(:,4:6) = Data;
Voronoi = VorArea3D(Data1);
S = 1./Voronoi{1};
%   Monte-Carlo, can be parallelized: for -> parfor, ...
tic;
iter = 50; %number of iterations
Nbins = round(1*(length(S))^(1/3)); %can be adjusted to smooth the histogram
lim = 2.7*median(S(S<Inf)); %can be adjusted to cover the needed range
[counts, centers] = hist(S(S<lim),Nbins);
binsize = centers(2) - centers(1);
counts_r = zeros(iter, Nbins);
for j = 1:iter
Rand = [rand(size(Data1,1),1)*1000, rand(size(Data1,1),1)*1000, rand(size(Data1,1),1)*1000];
ABrand = Data1;
ABrand(:,4:6) = Rand;

Voronoi_r = VorArea3D( ABrand );
Sr = 1./Voronoi_r{1};
[counts_r(j,:), ~] = hist(Sr(Sr<lim),centers);
end
toc;
for i = 1:2
    MeanCounts = mean(counts_r);
    StdCounts = std(counts_r);
end

%%   Plot curves, check the intersections manually and put them to the next section as 'thresh' and 'thresh2'
figure; plot(centers, counts, centers, MeanCounts + 1.96*StdCounts, '--m', centers, MeanCounts, '-.g', centers, MeanCounts - 1.96*StdCounts, '--m', 'LineWidth', 1.5); %1.96* to cover 95% of normal distribution, change if needed

set(gca, 'TickLength', [0.02,0.05], 'FontSize', 14);
legend('Data', 'Confidence envelope', 'Mean of randomized data');

%print('curvesr50bg4E-6cl40E-6.tif', '-dtiff', '-r300');
%   find intersection
%[xi, ~] = polyxpoly(centers, counts, centers, MeanCounts);
%%   Segment diagram
thresh = 6.769E3; % first intersection between the blue and the green lines
thresh2 = 2.073E4; %second intersection, put Inf if not needed
% %
% [ BW, clusters] = VoronoiSegmentation3D( Data1, thresh, 1, [0, Inf] );
% BW = BW(0:1200,0:1200);
% RGB = cat(3,BW,zeros(size(BW)),zeros(size(BW)));
% figure(); imshow(RGB);
% 
% imwrite(uint8(255*RGB), 'Segmented1x0.0025_0.000125.tif');

%
%figure; plot(centers, counts, centers, MeanCounts + 1.96*StdCounts, '--m', centers, MeanCounts, '-.g', centers, MeanCounts - 1.96*StdCounts, '--m', 'LineWidth', 1.5);
%Vav = (1000*1000*1000)/length(S);
%figure; plot(centers, counts * Vav / (length(S) * binsize), centers, (MeanCounts + 1.96*StdCounts) * Vav / (length(S) * binsize), '--m', centers, MeanCounts * Vav / (length(S) * binsize), '-.g', centers, (MeanCounts - 1.96*StdCounts) * Vav / (length(S) * binsize), '--m', centers, 3125/24*(centers/Vav).^4.*exp(-5*centers/Vav), '.k', 'LineWidth', 1.5);
figure; plot(centers, counts / (length(S) * binsize), '-b', centers, (MeanCounts + 1.96*StdCounts) / (length(S) * binsize), '--m', centers, MeanCounts / (length(S) * binsize), '-.g', centers, (MeanCounts - 1.96*StdCounts) / (length(S) * binsize), '--m', 'LineWidth', 2);
legend('Data', 'Confidence envelope', 'Mean of randomized distribution');
set(gca, 'TickLength', [0.01,0.01], 'FontSize', 16, 'LineWidth', 1);
xlabel('volume, nm^3','FontSize', 16);
ylabel('pdf','FontSize', 16);
title({'R_c_l=30nm; N=175(cl); d_c_l=4\times10^-^4; d_{bgd}=0.6\times10^-^4','Distibution of volumes of Voronoi polyhedrons'},'FontSize', 16);
xlim([0 3.975E4]);
ylim([0 6.5E-5]);
%% print curves, check the file name
print('curvesr30bg0.6E-4cl4E-4N175erf.tif', '-dtiff', '-r300');
%% keep points based on polyhedra
Acl = Voronoi{4}(S<thresh,:);
Colcl = repmat([0 0.5 0], size(Acl,1), 1);
Ard = Voronoi{4}(S<thresh2&S>thresh,:);
Colrd = repmat([0.85 0.85 0], size(Ard,1), 1);
Abg = Voronoi{4}(S>thresh2,:);
Colbg = repmat([0.85 0.85 0.85], size(Abg,1), 1);
Ashow = [Acl; Ard; Abg];
Colshow = [Colcl; Colrd; Colbg];
%% colorful scatter plot
figure; scatter3(Ashow(:,1), Ashow(:,2), Ashow(:,3),15, Colshow, '.'); axis equal;
view(20,10);
set(gca, 'TickLength', [0.01,0.01], 'FontSize', 14, 'LineWidth', 0.5, 'Xtick', 0:200:1000, 'Ytick', 0:500:1000, 'Ztick', 200:200:1000);
xlabel('nm','FontSize', 14);
ylabel('nm','FontSize', 14);
zlabel('nm','FontSize', 14);
%title({'R_c_l=30nm; N=175(cl); d_c_l=3\times10^-^4; d_{bgd}=0.6\times10^-^4','Distibution of volumes of Voronoi polyhedrons'},'FontSize', 16);
xlim([0 1000]);
ylim([0 1000]);
zlim([0 1000]);
%% print the scatter plot, check the filename
print('Col_r30bg0.6E-4cl4E-4N175erf.tif', '-dtiff', '-r400');
%% Interpolate to an image - can be slow!
tic;
[I, Avor] = drawVor3Dgeneral( Data1, 5, [1000,1000,1000]); %(data, voxelsize, field of view)
toc;
%% Binarize
InvThresh = 1/thresh;
BW = I;
BW = BW >= InvThresh;
%% Show the binarized volume
figure; patch(isosurface(I,InvThresh), 'FaceColor', [0.2 0.5 0.2],'FaceAlpha',.7,'LineStyle', 'none','FaceLighting', 'gouraud','BackFaceLighting','reverselit'); axis equal
view(20,10);
set(gca,'XGrid','on','YGrid','on','ZGrid','on', 'GridLineStyle',':','TickLength', [0.01,0.01], 'FontSize', 14, 'LineWidth', 0.5, 'Xtick', 0:40:200, 'Ytick', 0:100:200, 'Ztick', 40:40:200, 'XTickLabel', {'0','200','400','600','800','1000'}, 'YTickLabel', {'0','500','1000'}, 'ZTickLabel', {'200','400','600','800','1000'} );
xlabel('nm','FontSize', 14);
ylabel('nm','FontSize', 14);
zlabel('nm','FontSize', 14);
%title({'R_c_l=30nm; N=175(cl); d_c_l=3\times10^-^4; d_{bgd}=0.6\times10^-^4','Distibution of volumes of Voronoi polyhedrons'},'FontSize', 16);
xlim([0 200]);
ylim([0 200]);
zlim([0 200]);
light;
%% Print the binarized volume, check the filename
print('segm_r30bg0.6E-4cl4E-4N175erf.tif', '-dtiff', '-r400');
%% Find properties of the clusters
%imagesc(sum(BW(:,:,13:13),3));
stats = regionprops(BW, 'Area', 'Centroid', 'BoundingBox', 'SubarrayIdx', 'Image', 'FilledImage', 'FilledArea', 'PixelIdxList', 'PixelList');
% exclude Clusters smaller than a threshold
c = 0;
for i = 1:size(stats,1)
if stats(i-c).Area < 250 % Volume threshold
   stats(i-c) = [];
   c = c + 1;
end
end
Volumes = cat(1, stats.Area);
Volumes = 125* Volumes; %convert to nm^3 for voxelsize = 5nm
%% Draw equivalent radii
figure; hist(((3*Volumes)/(4*pi)).^(1/3),25);
xlim([3.1 34]);
ylim([0 340]);
xlabel('radius, nm','FontSize', 16);
ylabel('frequency','FontSize', 16);
title({'R_c_l=30nm; N=175(cl); d_c_l=4\times10^-^4; d_{bgd}=0.6\times10^-^4','Equivalent radius of clusters'},'FontSize', 16);
set(gca,'FontSize', 16, 'LineWidth',1);
%% save the plot, check filename
print('C:\Documents\MATLAB\Voronoi demostration\3D\Rr30bg0.6E-4cl4E-4N175erf.tif', '-dtiff', '-r300');
%% mean and standard deviation of the radii - show on screen
mean(((3*Volumes)/(4*pi)).^(1/3))
std(((3*Volumes)/(4*pi)).^(1/3))
%patch(isosurface(I,InvThresh))
%build a 3D histogram (size of BW) and count events from PixelIdxList
%% Nevents per cluster
%clear N
Ihist = drawHist3D(Data1, 5, [1000,1000,1000]);
for i = 1:size(stats,1)
PixelIdxList = stats(i).PixelIdxList;
N(i) = sum(Ihist(PixelIdxList));
end
mean(N)
std(N)