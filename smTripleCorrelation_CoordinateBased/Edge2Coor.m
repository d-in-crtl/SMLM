function [coordinates, coordinates_norm] = Edge2Coor(indx_13, indx_12, indx_23, seqmode, NormMode, normDist)

theta_1 = acos((indx_12.^2+indx_13.^2-indx_23.^2)./(2.*indx_12.*indx_13));
theta_2 = acos((indx_12.^2+indx_23.^2-indx_13.^2)./(2.*indx_12.*indx_23));
theta_3 = acos((indx_13.^2+indx_23.^2-indx_12.^2)./(2.*indx_13.*indx_23));

dise = 2 .* (round(rand(size(indx_13, 1), 1)) - 0.5);
        
switch seqmode
    case 123
        theta_1 = dise.*theta_1;
        
        x2 =  0.5 .* indx_12;
        y2 = x2 - x2;
        x1 = -0.5 .* indx_12;
        y1 = x2 - x2;
        x3 = indx_13 .* cos(theta_1) + x1;
        y3 = indx_13 .* sin(theta_1);
        
        switch NormMode
            case 'const'
                normfold = normDist ./ indx_12;
                x2_norm =  0.5 .* indx_12 .* normfold;
                y2_norm = x2_norm - x2_norm;
                x1_norm = -0.5 .* indx_12 .* normfold;
                y1_norm = x2_norm - x2_norm;
                x3_norm = indx_13 .* normfold .*cos(theta_1) + x1_norm;
                y3_norm = indx_13 .* normfold .*sin(theta_1);
            case 'perimeter'
                normfold = normDist ./ (indx_12 + indx_13 + indx_23);
                x2_norm =  0.5 .* indx_12 .* normfold;
                y2_norm = -0.5 .* indx_13 .* normfold .*sin(theta_1);
                x1_norm = -0.5 .* indx_12 .* normfold;
                y1_norm = -0.5 .* indx_13 .* normfold .*sin(theta_1);
                x3_norm = indx_13 .* normfold .*cos(theta_1) + x1_norm;
                y3_norm =  0.5 .* indx_13 .* normfold .*sin(theta_1);
            otherwise
                error('no normalization mode specified');
        end
        coordinates = [x1, y1, x2, y2, x3, y3]; 
        coordinates_norm = [x1_norm, y1_norm, x2_norm, y2_norm, x3_norm, y3_norm];
    
    case 132
        theta_1 = dise.*theta_1;
        
        x3 =  0.5 .* indx_13;
        y3 = x3 - x3;
        x1 = -0.5 .* indx_13;
        y1 = x3 - x3;
        x2 = indx_12 .* cos(theta_1) + x1;
        y2 = indx_12 .* sin(theta_1);
        
        switch NormMode
            case 'const'
                normfold = normDist ./ indx_13;
                x3_norm =  0.5 .* indx_13 .* normfold;
                y3_norm = x3_norm - x3_norm;
                x1_norm = -0.5 .* indx_13 .* normfold;
                y1_norm = x3_norm - x3_norm;
                x2_norm = indx_12 .* normfold .*cos(theta_1) + x1_norm;
                y2_norm = indx_12 .* normfold .*sin(theta_1);
            case 'perimeter'
                normfold = normDist ./ (indx_12 + indx_13 + indx_23);
                x3_norm =  0.5 .* indx_13 .* normfold;
                y3_norm = -0.5 .* indx_12 .* normfold .*sin(theta_1);
                x1_norm = -0.5 .* indx_13 .* normfold;
                y1_norm = -0.5 .* indx_12 .* normfold .*sin(theta_1);
                x2_norm = indx_12 .* normfold .*cos(theta_1) + x1_norm;
                y2_norm =  0.5 .* indx_12 .* normfold .*sin(theta_1);
            otherwise
                error('no normalization mode specified');
        end
        coordinates = [x1, y1, x3, y3, x2, y2]; 
        coordinates_norm = [x1_norm, y1_norm, x3_norm, y3_norm, x2_norm, y2_norm];
        
    case 213
        theta_2 = dise.*theta_2;
        
        x1 =  0.5 .* indx_12;
        y1 = x1 - x1;
        x2 = -0.5 .* indx_12;
        y2 = x1 - x1;
        x3 = indx_23 .* cos(theta_2) + x2;
        y3 = indx_23 .* sin(theta_2);
        
        switch NormMode
            case 'const'
                normfold = normDist ./ indx_12;
                x1_norm =  0.5 .* indx_12 .* normfold;
                y1_norm = x1_norm - x1_norm;
                x2_norm = -0.5 .* indx_12 .* normfold;
                y2_norm = x1_norm - x1_norm;
                x3_norm = indx_23 .* normfold .*cos(theta_2) + x2_norm;
                y3_norm = indx_23 .* normfold .*sin(theta_2);
            case 'perimeter'
                normfold = normDist ./ (indx_12 + indx_13 + indx_23);
                x1_norm =  0.5 .* indx_12 .* normfold;
                y1_norm = -0.5 .* indx_23 .* normfold .*sin(theta_2);
                x2_norm = -0.5 .* indx_12 .* normfold;
                y2_norm = -0.5 .* indx_23 .* normfold .*sin(theta_2);
                x3_norm = indx_23 .* normfold .*cos(theta_2) + x2_norm;
                y3_norm =  0.5 .* indx_23 .* normfold .*sin(theta_2);
            otherwise
                error('no normalization mode specified');
        end
        coordinates = [x2, y2, x1, y1, x3, y3]; 
        coordinates_norm = [x2_norm, y2_norm, x1_norm, y1_norm, x3_norm, y3_norm];
        
    case 231
        theta_2 = dise.*theta_2;
        
        x3 =  0.5 .* indx_23;
        y3 = x3 - x3;
        x2 = -0.5 .* indx_23;
        y2 = x3 - x3;
        x1 = indx_12 .* cos(theta_2) + x2;
        y1 = indx_12 .* sin(theta_2);
        
        switch NormMode
            case 'const'
                normfold = normDist ./ indx_23;
                x3_norm =  0.5 .* indx_23 .*normfold;
                y3_norm = x3_norm - x3_norm;
                x2_norm = -0.5 .* indx_23 .*normfold;
                y2_norm = x3_norm - x3_norm;
                x1_norm = indx_12 .*normfold .*cos(theta_2) + x2_norm;
                y1_norm = indx_12 .*normfold .*sin(theta_2);
            case 'perimeter'
                normfold = normDist ./ (indx_12 + indx_13 + indx_23);
                x3_norm =  0.5 .* indx_23 .*normfold;
                y3_norm = -0.5 .* indx_12 .*normfold .*sin(theta_2);
                x2_norm = -0.5 .* indx_23 .*normfold;
                y2_norm = -0.5 .* indx_12 .*normfold .*sin(theta_2);
                x1_norm = indx_12 .*normfold .*cos(theta_2) + x2_norm;
                y1_norm =  0.5 .* indx_12 .*normfold .*sin(theta_2);
            otherwise
                error('no normalization mode specified');
        end
        coordinates = [x2, y2, x3, y3, x1, y1]; 
        coordinates_norm = [x2_norm, y2_norm, x3_norm, y3_norm, x1_norm, y1_norm];
        
    case 312
        theta_3 = dise.*theta_3;
        
        x1 =  0.5 .* indx_13;
        y1 = x1 - x1;
        x3 = -0.5 .* indx_13;
        y3 = x1 - x1;
        x2 = indx_23 .* cos(theta_3) + x3;
        y2 = indx_23 .* sin(theta_3);
        
        switch NormMode
            case 'const'
                normfold = normDist ./ indx_13;
                x1_norm =  0.5 .* indx_13 .* normfold;
                y1_norm = x1_norm - x1_norm;
                x3_norm = -0.5 .* indx_13 .* normfold;
                y3_norm = x1_norm - x1_norm;
                x2_norm = indx_23 .* normfold .*cos(theta_3) + x3_norm;
                y2_norm = indx_23 .* normfold .*sin(theta_3);
            case 'perimeter'
                normfold = normDist ./ (indx_13 + indx_12 + indx_23);
                x1_norm =  0.5 .* indx_13 .* normfold;
                y1_norm = -0.5 .* indx_23 .* normfold .*sin(theta_3);
                x3_norm = -0.5 .* indx_13 .* normfold;
                y3_norm = -0.5 .* indx_23 .* normfold .*sin(theta_3);
                x2_norm = indx_23 .* normfold .*cos(theta_3) + x3_norm;
                y2_norm = 0.5 .* indx_23 .* normfold .*sin(theta_3);
            otherwise
                error('no normalization mode specified');
        end
        coordinates = [x3, y3, x1, y1, x2, y2]; 
        coordinates_norm = [x3_norm, y3_norm, x1_norm, y1_norm, x2_norm, y2_norm];
        
    case 321
        theta_3 = dise.*theta_3;
        
        x2 =  0.5 .* indx_23;
        y2 = x2 - x2;
        x3 = -0.5 .* indx_23;
        y3 = x2 - x2;
        x1 = indx_13 .* cos(theta_3) + x3;
        y1 = indx_13 .* sin(theta_3);
        
        switch NormMode
            case 'const'
                normfold = normDist ./ indx_23;
                x2_norm =  0.5 .* indx_23 .* normfold;
                y2_norm = x2_norm - x2_norm;
                x3_norm = -0.5 .* indx_23 .* normfold;
                y3_norm = x2_norm - x2_norm;
                x1_norm = indx_13 .* normfold .* cos(theta_3) + x3_norm;
                y1_norm = indx_13 .* normfold .* sin(theta_3);
            case 'perimeter'
                normfold = normDist ./ (indx_23 + indx_13 + indx_12);
                x2_norm =  0.5 .* indx_23 .* normfold;
                y2_norm = -0.5 .* indx_13 .* normfold .* sin(theta_3);
                x3_norm = -0.5 .* indx_23 .* normfold;
                y3_norm = -0.5 .* indx_13 .* normfold .* sin(theta_3);
                x1_norm = indx_13 .* normfold .* cos(theta_3) + x3_norm;
                y1_norm =  0.5 .* indx_13 .* normfold .* sin(theta_3);
            otherwise
                error('no normalization mode specified');
        end
        coordinates = [x3, y3, x2, y2, x1, y1]; 
        coordinates_norm = [x3_norm, y3_norm, x2_norm, y2_norm, x1_norm, y1_norm];
    
    otherwise
        error('colors do not match');
end

