% This script takes global, 240 year fields of temperature and carbon
% (available upon request), and selects a section which approximates the
% location of the RAPID array. It will then use 100 years of control run
% data in order to generate a kappa_r field, as well as loading this
% section directly from computed values.

load('/home/ct/MATLAB/BIGMATFILES/DIC_RAD.mat');
load('/home/ct/MATLAB/BIGMATFILES/DIC_CTR.mat');

DIC_CTR = single(DIC_CTR); DIC_RAD = single(DIC_RAD);

dCnat = DIC_RAD - DIC_CTR;


% Now we need to get dCanth too.

load('/home/ct/MATLAB/BIGMATFILES/DIC_COU.mat');

DIC_COU = single(DIC_COU);

dCanth = DIC_COU - DIC_RAD;

[dCnat, dCanth, DIC_COU, DIC_CTR, DIC_RAD] = globalToRapid(dCnat, dCanth, DIC_COU, DIC_CTR, DIC_RAD);

load('/home/ct/MATLAB/BIGMATFILES/TMP_CTR.mat');
load('/home/ct/MATLAB/BIGMATFILES/TMP_COU.mat');

TMP_CTR = single(TMP_CTR); TMP_COU = single(TMP_COU);
[TMP_CTR,TMP_COU] = globalToRapid(TMP_CTR,TMP_COU);

save RAPID_TMP_DIC_Fields.mat dCnat dCanth DIC_COU DIC_CTR DIC_RAD TMP_CTR TMP_COU
%%
load RAPID_TMP_DIC_Fields.mat


% Now load kappa_r values in. This .mat file also contains polynomial fits
% which were used as a sanity check.
% This file is generated in TestingShortAndLongTermGlobalOcean.m
load('./KappaRandPolyfit.mat','kappa_s','ecc');
kappa_r = kappa_s .* ecc;
kappa_r = globalToRapid(kappa_r);

%

kappa_s_Calculated = NaN(76,64);
ecc_Calculated = NaN(76,64);
p1Vals_Calculated = NaN(76,64);


for i = 1:76
    for j = 1:64
            x = squeeze(DIC_CTR(i,j,1:240));
            y = squeeze(TMP_CTR(i,j,1:240));
            [kappa_s_Calculated(i,j),ecc_Calculated(i,j)] = compute_kappa_r(x,y,0,10,'ecc');
            p1 = polyfit(x,y,1);
            p1Vals_Calculated(i,j) = p1(1);
    end
    fprintf('Done on longitude box %d\n',i);
end
kappa_r_Calculated = kappa_s_Calculated .* ecc_Calculated;
% Show that the calculated and saved versions are the same

pcolor(kappa_r - kappa_r_Calculated); shading flat
%% Now we can calculate excess and redistributed temperature from 
% these fields.

load gammaValues.mat

smoothedGammaValuesArr = repmat(smoothedGammaValues,[1 76 64]);
smoothedGammaValuesArr = permute(smoothedGammaValuesArr,[2 3 1]);

AdjCnat = dCnat + smoothedGammaValuesArr .* dCanth;

dTMPr = kappa_r .* AdjCnat;
dTMP = TMP_COU - TMP_CTR;
dTMPe = dTMP - dTMPr;

% And plot this off
figure; 
subplot(1,3,1);
pcolor(nanmean(dTMP(:,:,231:240),3)');
set(gca,'Ydir','reverse'); shading flat;
title('Total Temperature Change');

subplot(1,3,2);
pcolor(nanmean(dTMPr(:,:,231:240),3)');
set(gca,'Ydir','reverse'); shading flat;
title('Redistributed Temperature');

subplot(1,3,3);
pcolor(nanmean(dTMPe(:,:,231:240),3)');
set(gca,'Ydir','reverse'); shading flat;
title('Excess Temperature');

