function [outputGammaValue] = calculateGammaValue(CnatField,CanthField,kappa_r,CarbonOrTemp)
% This function will take a global (3D or 4D if we include multiple times)
% Cnat Field, Canth Field, a global kappa_r field, and finally a switch 
% (CarbonOrTemp), which allows you to specify which quantity we want to 
% enforce a globally integrated value of zero. Set CarbonOrTemp = 1 to
% enforce globally integrated redistributed heat == 0, 0 for globally
% integrated redistributed carbon == 0.


if size(CnatField) ~= size(CanthField)
    error('Canth and Cnat fields must be the same size');
end

try
    dV = evalin('base','dV');
catch
    load grid_data_areas_masks.mat dV
end

outputGammaValue = NaN(1);

if ndims(CnatField) ~= ndims(kappa_r)
    nDimField = ndims(kappa_r);
    truthVector = NaN(nDimField,1);
    for i = 1:nDimField
      truthVector(i) = (size(CnatField,i) == size(kappa_r,i));
    end
    
    if nanmean(truthVector ~= 1)
        error('If carbon fields have a time dimension, it must be the last dimension. If time dimension is the last dimension, check spatial dimensions agree');
    end
    outputGammaValue = NaN(size(CnatField,nDimField+1),1);
end

% The search is run as follows: we will look for an optimal gamma value to
% 3 decimal places. In order to speed this search up, we will first look
% for minima to one decimal place. A second search follows to two decimal
% places, followed by a 3rd search to 3 decimal places. 
% For example: Search 1 yields a minimum value of 0.3, searching the range
% 0:0.1:1. Search 2 will then search the range 0.2:0.01:0.4, yielding an
% optimal value of 0.24. Search 3 will then search the range
% 0.23:0.001:0.24. 
% This search technique assumes that there is a unique minimum globally 
% integrated redistributed temperature or carbon value as a function of
% gamma. As Canth ought to be positive definite (excluding potential
% machine accuracy issues), this should be the case, or at least an
% extremely good approximation.


% Single precision to circumvent memory limit.
CnatField = single(CnatField); 
CanthField = single(CanthField);

% kappa_r is just there to convert from carbon change to temperature change
% If we are minimising carbon, we can just set it to one to avoid rewriting
% code
if CarbonOrTemp == 0
    kappa_r = ones(size(kappa_r));
end

parfor i = 1:length(outputGammaValue) % Need to repeat for each time point
    fprintf('Calculating gamma value for year %d of run\n',i);
    valRange = 0:0.1:1;
    Residuals = NaN(size(valRange));
    CnatFieldOneYear = CnatField(:,:,:,i);
    CanthFieldOneYear= CanthField(:,:,:,i); 
    for j = 1:length(valRange)
        Residuals(j) = abs( nansum(kappa_r(:) .* dV (:) .* ( CnatFieldOneYear(:) + valRange(j) * CanthFieldOneYear(:)) ));
    end

    [~,MinValue] = min(Residuals);

    valRange = valRange(MinValue)-0.1:0.01:valRange(MinValue)+0.1;
    valRange(valRange<0 | valRange>1) = [];
    
    for j = 1:length(valRange)
        Residuals(j) = abs( nansum(kappa_r(:) .* dV (:) .*  ( CnatFieldOneYear(:) + valRange(j) * CanthFieldOneYear(:)) ));
    end

    [~,MinValue] = min(Residuals);
        
    valRange = valRange(MinValue)-0.01:0.001:valRange(MinValue)+0.01;
    valRange(valRange<0 | valRange>1) = [];
    
    for j = 1:length(valRange)
        Residuals(j) = abs( nansum(kappa_r(:) .* dV (:) .*  ( CnatFieldOneYear(:) + valRange(j) * CanthFieldOneYear(:)) ));
    end

    [~,MinValue] = min(Residuals);
    
    outputGammaValue(i) = valRange(MinValue);
end
    

