function [triple_trans, results, coordinates, coordinates_norm] = TripleRead(triplemat, PC_1, PC_2, PC_3, boxradii, dim, rho_res, NormMode, normDist, TripleCalSeq, TripleDispSeq)
% ========================= Sequence Info =================================
seq = zeros(1,3);
seq(1) = find(TripleCalSeq == TripleDispSeq(1));
seq(2) = find(TripleCalSeq == TripleDispSeq(2));
seq(3) = find(TripleCalSeq == TripleDispSeq(3));
seqmode = seq(1)*100 + seq(2)*10 + seq(3);

% ========================== Theta Info ===================================
triplemat = permute(triplemat, [1 3 2]);
triplemat(dim+1:end,:,:) = [];
triplemat(:,dim+1:end,:) = [];

SizeofTriple = size(triplemat);
if dim > min(SizeofTriple(1), SizeofTriple(2))
    error('dim is bigger than the data cube itself');
end

rad_res = 2*pi/SizeofTriple(3);
rad_center = ((1:SizeofTriple(3))' - 0.5).*rad_res - pi;

% ======================== Extracting CrossCorrelation ====================

for ii = 1 : SizeofTriple(1)
    edge_1 = PC_1(ii, 1);
    for jj = 1 : SizeofTriple(2)
        edge_2 = PC_2(jj, 1);
        for kk = 1 : SizeofTriple(3)
            theta = rad_center(kk);
            edge_3 = sqrt(edge_1^2 + edge_2^2 - 2*edge_1*edge_2*cos(theta));
            trans_kk = ceil(edge_3 / rho_res);
            triplemat(ii, jj, kk) = triplemat(ii, jj, kk) - (PC_1(ii, 2)-1) - (PC_2(jj, 2)-1) - (PC_3(trans_kk, 2)-1);
        end
    end
end
       
% ==================== building up transforming cube ======================
triple_trans = ones(dim, dim, dim);

trans_cube_1 = repmat((1:dim)', [1, dim, SizeofTriple(3)]);
trans_cube_2 = repmat((1:dim),  [dim, 1, SizeofTriple(3)]);
trans_cube_3 = repmat(permute(rad_center, [2 3 1]), [dim, dim, 1]);
trans_cube = sqrt(trans_cube_1.^2 + trans_cube_2.^2 - 2.*trans_cube_1.*trans_cube_2.*cos(trans_cube_3));

% ====================== transforming and smoothing ======================= 
rho_edge = (0 : dim)'; 
for d1 = 1 : dim
    for d2 = 1 : dim
        tmp_vt = permute(trans_cube(d1, d2, :), [3 2 1]);
        [~, ~, bin] = histcounts(tmp_vt, rho_edge);
        for d3 = 1 : dim
            if ~isempty(triplemat(d1, d2, bin == d3))
                triple_trans(d1, d2, d3) = mean(triplemat(d1, d2, bin == d3));
            end
        end
    end
end

triple_trans = imgaussfilt3(triple_trans, boxradii);

evaluate = triple_trans(:);
evaluate(evaluate <= 1) = [];

threshold = mean(evaluate) + 2.5 * std(evaluate);
disp(['mean = ' num2str(mean(evaluate)) ', std = ' num2str(std(evaluate)) ', max = ' num2str(max(evaluate)) ', thresh = ' num2str(threshold)]);
%triple_trans(triple_trans <= threshold) = 1;

% ================== find the local maximum of the cube ===================
local_mask = imregionalmax(triple_trans) & triple_trans > threshold;
indx = find(local_mask);

I = triple_trans(indx);
        
indx_23 = floor((indx - 1)./(dim*dim)) + 1; % dist 2_3
indx_12_13 = mod((indx - 1), (dim*dim)) + 1; 
indx_12 = floor((indx_12_13 - 1)./dim) + 1; % dist 1_2
indx_13 = mod((indx_12_13 - 1), dim) + 1; % dist 1_3

% ==================== filter the rough local maximums ====================
mask = zeros(size(indx, 1), 1);
for k = 1 : size(indx, 1)
    subvolume = triple_trans(max(indx_13(k)-boxradii(1),1):min(indx_13(k)+boxradii(1),dim), max(indx_12(k)-boxradii(2),1):min(indx_12(k)+boxradii(1),dim), max(indx_23(k)-boxradii(3),1):min(indx_23(k)+boxradii(3),dim));
    if triple_trans(indx_13(k), indx_12(k), indx_23(k)) < max(subvolume(:)) || max(subvolume(:)) <= 1
        mask(k) = 1;
        continue
    end
    corner_1 = max(indx_13(k) - boxradii(1), 1); % corner 1_3
    corner_2 = max(indx_12(k) - boxradii(2), 1); % corner 1_2
    corner_3 = max(indx_23(k) - boxradii(3), 1); % corner 2_3
            
    [volsz1, volsz2, volsz3] = size(subvolume);
    [subcoor_2, subcoor_1, subcoor_3] = meshgrid((0.5:volsz2-0.5), (0.5:volsz1-0.5), (0.5:volsz3-0.5));
    subcenter_1 = sum(sum(sum(subvolume.*subcoor_1)))/sum(subvolume(:));
    subcenter_2 = sum(sum(sum(subvolume.*subcoor_2)))/sum(subvolume(:));
    subcenter_3 = sum(sum(sum(subvolume.*subcoor_3)))/sum(subvolume(:));
            
    indx_13(k) = corner_1-1 + subcenter_1;
    indx_12(k) = corner_2-1 + subcenter_2;
    indx_23(k) = corner_3-1 + subcenter_3;
            
    if round(indx_13(k)) < 0 || round(indx_13(k)) > dim || round(indx_12(k)) < 0 || round(indx_12(k)) > dim || round(indx_23(k)) < 0 || round(indx_23(k)) > dim
        mask(k) = 1;
        continue
    end
    
    %if round(indx_13(k)) < 15 && round(indx_12(k)) < 15 && round(indx_23(k)) < 15
    %    mask(k) = 1;
    %    continue
    %end
   
end
        
% =============================== rejections ==============================
rej1 = mask > 0;

indx_23(rej1) = []; % dist 2_3
indx_12(rej1) = []; % dist 1_2
indx_13(rej1) = []; % dist 1_3
I(rej1) = [];

cos_theta_1 = (indx_12.^2+indx_13.^2-indx_23.^2)./(2.*indx_12.*indx_13);
cos_theta_2 = (indx_12.^2+indx_23.^2-indx_13.^2)./(2.*indx_12.*indx_23);
cos_theta_3 = (indx_13.^2+indx_23.^2-indx_12.^2)./(2.*indx_13.*indx_23);

rej2 = abs(cos_theta_1) > 1 | abs(cos_theta_2) > 1 | abs(cos_theta_3) > 1;

indx_23(rej2) = []; % dist 2_3
indx_12(rej2) = []; % dist 1_2
indx_13(rej2) = []; % dist 1_3
I(rej2) = [];

results = [indx_13.*rho_res, indx_12.*rho_res, indx_23.*rho_res, I];

% ======================== calculate the coordinates ======================

[coordinates, coordinates_norm] = Edge2Coor(indx_13.*rho_res, indx_12.*rho_res, indx_23.*rho_res, seqmode, NormMode, normDist);

