function [ Image ] = drawHist3D( A, p, fov, ko, ke )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here.
% in nm
if isstruct(A)
A = A.data;
end
if ~exist('ko', 'var')
    ko = 1;
end
if ~exist('ke', 'var')
    ke = size(A, 1);
end
if ~exist('fov', 'var') | fov == 0
    fov = [max(A(:,4)) max(A(:,5)) max(A(:,6))];
end

A = A(ko:ke, :);

NZ = round(fov(3)/p); % number of slices in Z
Image = zeros(floor(fov(1)/p), floor(fov(2)/p), NZ);
l = size(A,1); % number of rows
n = ones(NZ,1);
Anew = cell(NZ,1);
for k = 1:l
    for i = 1:NZ;
    if (A(k,6) > p * (i - 1)) && (A(k,6) <= p * i)
        Anew{i}(n(i), :) = A(k,:);
        n(i) = n(i) + 1;
    end
    end
end
%parfor i = 1:NZ
for i = 1:NZ
    if ~isempty(Anew{i})
        I = drawHist2D(Anew{i}, p, fov(1:2));
        Image(:,:,i) = I;
    else
        Image(:,:,i) = zeros(floor(fov(1:2)/p));
    end
end
end