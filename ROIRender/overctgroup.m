function [out] = overctgroup(in, var, frm)
% in = [x_in, y_in]; matrix unit: px 
% precision = var of the corresponding x_in and y_in; vector unit: px^2
% frm = frame number of corresponding x_in and y_in

% out = [x_out, y_out, var_out]; unit: px

last_frm = max(frm(:));
first_frm = min(frm(:));

out = [0, 0, 0];

modified_x = in(frm == first_frm, 1);
modified_y = in(frm == first_frm, 2);
modified_var = var(frm == first_frm, 1);

for frm_idx = first_frm : last_frm - 1
    
    current_x = modified_x;
    current_y = modified_y;
    current_var = modified_var;
    
    next_x = in(frm == frm_idx+1, 1);
    next_y = in(frm == frm_idx+1, 2);
    next_var = var(frm == frm_idx+1, 1);
    
    var_tmp = cat(1, current_var, next_var);
    sigma = min(mean(sqrt(var_tmp(:, 1)), 1), median(sqrt(var_tmp(:, 1)), 1));
        
    dist = pdist2([next_x, next_y], [current_x, current_y]);
    
    mask = dist < 5*sigma;
    
    out_idx = find(sum(mask, 1)' == 0);
    if ~isempty(out_idx)
        x_out = current_x(out_idx);
        y_out = current_y(out_idx);
        var_out = current_var(out_idx);
        out = cat(1, out, [x_out, y_out, var_out]);
    end
    
    modified_x = next_x;
    modified_y = next_y;
    modified_var = next_var;
    
    modify_idx = find(sum(mask, 2) > 0);
    if ~isempty(modify_idx)
        for i = 1 : size(modify_idx)
            tmp_idx = find(mask(modify_idx(i), :) > 0);
            tmp_idx = tmp_idx';
            modified_x(modify_idx(i)) = ( next_x(modify_idx(i))/next_var(modify_idx(i))+sum(current_x(tmp_idx)./current_var(tmp_idx), 1) ) ./ (1/next_var(modify_idx(i)) + sum(1./ current_var(tmp_idx), 1));
            modified_y(modify_idx(i)) = (next_y(modify_idx(i))/next_var(modify_idx(i)) + sum(current_y(tmp_idx)./current_var(tmp_idx), 1)) ./ (1/next_var(modify_idx(i)) + sum(1./ current_var(tmp_idx), 1));
            modified_var(modify_idx(i)) = (1+size(tmp_idx, 1))/(1/next_var(modify_idx(i)) + sum(1./current_var(tmp_idx), 1));
        end
    end
end
out(1,:) = [];
