%% Function to generate affinity matrix and distinctness constraint matrix
%
% Inputs:
%       - A:            Set of associations
%       - P1,P2:        Point clouds
%
% Output:
%       - M:            Affinity matrix
%       - C:            Constraint matrix (distinctness & consistency)
%
%
function [M, C] = affinity_matrix(A, P1,P2, eps, sig)

nA = size(A,1); % Total number of associations

% Preallocate matrices (symmetric, use sparse for speed?)
M = eye(nA)*0.5; % M will be added to its transpose, so diagonals will become 1
C = tril(true(nA)) - 0.5*eye(nA);  % C will be added to its transpose, so diagonal will become 1

% list all combinations of associations (in vector)
disp('Generating combination of associations...');
comb = zeros(nA*(nA-1)/2, 2);
itr = 0;
for i = 1 : nA-1
    jcol = (i+1:nA).';
    icol = i*ones(nA-i,1);
    comb(itr+1 : itr+(nA-i), :) = [jcol, icol];
    itr = itr+(nA-i);
%     for j = i+1 : nA
%         itr = itr + 1;
%         comb(itr,:) = [j,i];
%     end
end

disp('Computing consistency scores...');
% end points of the line segments corresponding to association pairs
p1_a = P1(A(comb(:,1), 1), :);
p1_b = P1(A(comb(:,2), 1), :);

% end points of the line segments corresponding to association pairs
p2_a = P2(A(comb(:,1), 2), :);
p2_b = P2(A(comb(:,2), 2), :);


% length of line segments
d1 = sqrt( sum( (p1_a - p1_b).^2 ,2) );
d2 = sqrt( sum( (p2_a - p2_b).^2 ,2) );

% residual distance
resid = abs(d1-d2);

Mvec = zeros(size(resid)); % lower triangular part of M matrix
Cvec = ones(size(resid)); % lower triangular part of C matrix

idx = (resid < eps); % index of consistent associations
Mvec(idx) = exp(-resid(idx).^2. / (2 * sig^2)); % consistency score
Cvec(~idx) = 0; % inconsistent associations

% inconsistent associations that start/end at same points (not one-to-one)
idx_z1 = (d1 == 0);
idx_z2 = (d2 == 0);
Mvec(idx_z1) = 0;
Mvec(idx_z2) = 0;
Cvec(idx_z1) = 0;
Cvec(idx_z2) = 0;

% put vectorized values into the lower triangular part
M( tril( (true(nA)-eye(nA))>0.5 ) ) = Mvec;
C( tril( (true(nA)-eye(nA))>0.5 ) ) = Cvec;

% Make matrices symmetric
M = M + M.';
C = C + C.';

disp('Affinity matrix constructed.');


