% This script will generate a kappa_r value for each point in the ORCA1
% grid as well as a polynomial fit value at each point which maybe be used
% as a sanity check


load('/home/ct/MATLAB/BIGMATFILES/DIC_CTR.mat');
load('/home/ct/MATLAB/BIGMATFILES/TMP_CTR.mat');

kappa_r = NaN(362,292,64);
ecc = NaN(362,292,64);
p1Vals = NaN(362,292,64);
%
for i = 1:362
    for j = 1:292
        for k = 1:64
            x = squeeze(DIC_CTR(i,j,k,1:100));
            y = squeeze(TMP_CTR(i,j,k,1:100));
            [kappa_r(i,j,k),ecc(i,j,k)] = compute_kappa_r(x,y,1,10,'ecc');
            p1 = polyfit(x,y,1);
            p1Vals(i,j,k) = p1(1);
        end
    end
    fprintf('Done on box %d\n',i);
end

%save KappaRandPolyfitValues.mat kappa_r ecc p1Vals
