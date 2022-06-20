function [I] = drawHist2D(A, p, fov)
% in nm

if ~exist('fov', 'var') | fov == 0
    fov = [max(A(:,4)) max(A(:,5))];
end
Anew(:,1) = A(:,5);
Anew(:,2) = A(:,4);
edges = {0:p:fov(1)-p, 0:p:fov(2)-p};
I = hist3 (Anew, 'Edges', edges);
end

