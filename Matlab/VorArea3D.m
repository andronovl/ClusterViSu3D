function [ Voronoi ] = VorArea3D( A1, ko, ke )
% Performs Voronoi tessellation on super-res data and calculates the 
% inverted areas for Voronoi cells
% Optional parameters: ko, ke

if isstruct(A1)
A1 = A1.data;
end
if ~exist('ko', 'var')
    ko = 1;
end
if ~exist('ke', 'var')
    ke = size(A1, 1);
end

Anew = A1(ko:ke, 4:6);
[A, ~, ic] = unique(Anew,'rows','stable');
% dt = delaunayTriangulation(A);
% [V, C] = voronoiDiagram(dt);
[V,C] = voronoin(A);
v = zeros(size(C,1),1);
for i = 1:size(C,1)
X = [V(C{i}, 1), V(C{i}, 2), V(C{i}, 3)];
if ~any(any(X == Inf))
    [~, v(i)] = convhulln(X, {'Qt','Qs'}); %volumes of polyhedra 
else
    v(i) = Inf;
end    
end
v = 1./v; %local density (nm^-3)
[n, ~] = histc(ic, unique(ic)); % n = number of repeats for each point with index ic.
v = v .* n;
v(isnan(v) | isinf(v)) = 0;
Voronoi = { v, V, C, A, ic };
