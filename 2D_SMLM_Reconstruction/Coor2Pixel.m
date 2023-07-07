function [out] = Coor2Pixel(results, CameraPxsz, targetsz, imszx, imszy, option)
% INPUT:  results should be the 'reports' given by the SMLM_sCMOS.m function, which inlcuds all the reconstruction information
%         results = [result_x, result_y, tot_sigma, result_precision, tot_photon, tot_bg, tot_pval, frm]
%         zm is the zooming factor, given by camera px size / target px size
%         imszy, imszx: original image height and width
%         options: 'b': reconstruct to binary image
%                  'I': reconstruct with fitted photon num
%                  's': reconstruct with fitted sigma
%                  'f': reconstruct with frame numbers
% OUTPUT: out: reconstructed image

x_edge = (0 : ceil(imszx * CameraPxsz / targetsz))' .* targetsz;
y_edge = (0 : ceil(imszy * CameraPxsz / targetsz))' .* targetsz;

x = results(:, 1) .* CameraPxsz;
y = results(:, 2) .* CameraPxsz;
photon = results(:, end-6);
frm = results(:, end-3);
llr = results(:, end);

[N, ~, ~, y_bin, x_bin] = histcounts2(y, x, y_edge, x_edge);
ind = N > 0;

if ~isempty(x)
    
    out = zeros(size(N));
    if option == 'b'
        out(ind) = 1;
    else
        if option == 'n'
            out = N;
        else
            for i = 1 : size(x, 1)
                if y_bin(i) ~= 0 && x_bin(i) ~= 0
                    switch option
                        case 'I'
                            out(y_bin(i), x_bin(i)) = out(y_bin(i), x_bin(i)) + photon(i);
                        case 'z'
                            out(y_bin(i), x_bin(i)) = out(y_bin(i), x_bin(i)) + sigma(i);
                        case 'llr'
                            out(y_bin(i), x_bin(i)) = out(y_bin(i), x_bin(i)) + llr(i);
                        case 'f'
                            out(y_bin(i), x_bin(i)) = frm(i); % put the last frame if this rendered pixel has more than 1 spots. 
                        otherwise
                            error('please choose an reconstruction option: b: binary, I: intensity, s: sigma, f: frame number ');
                    end
                end
            end
            
            out(ind) = out(ind) ./ N(ind);
            if option == 'z'
                out(ind) = out(ind) - min(out(:)) + 100;
            end
        end
    end
end
out = uint16(out);