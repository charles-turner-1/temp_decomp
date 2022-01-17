function [output,varargout] = compute_kappa_r(Cnat_timeseries,tmp_timeseries,verbosity,periodicity,varargin)

% Cnat_timeseries and tmp_timeseries must be vectors of the same length. This function will only work 
% out kappa_r at a single point for each call.
% verbosity must be either 1 or 0. Setting verbosity = 1 will give you a lot of feedback about what 
% the function is doing, setting 0 will result in no feedback. This is handy when running a lot of 
% points.
% periodicity can be set to whatever you like, but the standard usage I've been following is to set
% it to 12 if using monthly model data or 10 if using decadal data. It will determine how to block up 
% the timeseries to remove a mean from each value. Ie if you set periodicity = 12 with monthly data, 
% the timeseries will be blocked into years, and the yearly mean of each year removed to give monthly 
% - yearly mean data.
% varargin: must be 'det', 'ecc', or 'both'. If not specified, it will default to 'det'. This allows
% you to specify what kind of scaling factor is applied to account for fit uncertainty: determinant 
% or eccentricity. These are given as separate outputs.

n_vargin = numel(varargin);
n_vargout = nargout - 1;

if n_vargin & ~n_vargout
    warning('You have specified to return kappa_r and a scale factor, but function is only returning kappa_r. This is likely to cause confusion');
end

if n_vargin & n_vargin > 1
    error('Only 3 to 5 inputs can be entered');
elseif n_vargin == 0 & verbosity
    disp('No scale factor specified. Defaulting to determinant');
end

if isempty(periodicity)
    if verbosity
        disp('No periodicity specified. Assuming you want to use monthly data.')
        disp('Setting periodicity = 12');
    end
    periodicity = 12;
end

% This next condition should (touch wood) build in backwards compatability
% with previous fixed periodicity versions of the code.

if ~isnumeric(periodicity) | length(periodicity) ~=1
    if verbosity
        disp('Periodicity is not a number. Assuming backwards compatability cockup');
        disp('Assuming you want to use monthly data and setting periodicity = 12');
    end
    periodicity = 12;
end



if n_vargin == 1
    if  ~strcmpi(varargin{1}, 'ecc') 
        disp('If varargin specified, scale factor must be eccentricity');
        disp("Specify varargin == 'ecc' for eccentricity output")
        error('Scale factor specified incorrectly');
    end
end


Cnat_timeseries = squeeze(Cnat_timeseries); % Make sure we don't have extraneous dimensions 
tmp_timeseries = squeeze(tmp_timeseries); % Make sure we don't have extraneous dimensions

if ~isvector(tmp_timeseries) | ~isvector(Cnat_timeseries)
    error('Carbon, Temperature data must be vectors. Error in monthly data');
end

if size(Cnat_timeseries,2) ~= 1
    Cnat_timeseries = Cnat_timeseries'; % Make sure input data is a column vector
end

if size(tmp_timeseries,2) ~= 1
    tmp_timeseries = tmp_timeseries'; % Make sure input data is a column vector
end


if size(tmp_timeseries) ~= size(Cnat_timeseries)
    error('Carbon, Temperature vectors must be the same size. Error in monthly data.')
end

if ~any(~isnan(tmp_timeseries)) & ~any(~isnan(Cnat_timeseries))
    if verbosity
        disp('Temperature and Cnat vectors are all NaN. Assigning NaN output');
    end
    output = NaN;
    for i = 1:n_vargout
        varargout{i} = NaN;
    end
    return
end

if nanmean(Cnat_timeseries < 50) & nanmean(tmp_timeseries > 1500)
   [Cnat_timeseries,tmp_timeseries] = deal(tmp_timeseries,Cnat_timeseries);
   if verbosity
       disp("Tmp, Cnat seem to be misplaced. Switching")
   end
end

nan_idx = find(isnan(tmp_timeseries) | isnan(Cnat_timeseries));
if nan_idx
    tmp_timeseries(nan_idx) = 0;
    Cnat_timeseries(nan_idx) = 0;
end


high_frequency_tmp = compute_high_frequencies(tmp_timeseries,periodicity,verbosity);
high_frequency_Cnat = compute_high_frequencies(Cnat_timeseries,periodicity,verbosity);



% In some cases (ie. drift correction), we expect to see a period (ie.
% 30yrs) of temperature and Cnat being zero. Removing these has no effect
% as I've used PCA, but chuck out an warning anyway.

nanlist = find(high_frequency_Cnat==0 & high_frequency_tmp==0);
high_frequency_Cnat(nanlist) = [];
high_frequency_tmp(nanlist) = [];

if verbosity
    l = length(nanlist);
    if l
        warning(['Input data truncated by ',num2str(l),' data points due to zeros in input.'])
    end
    % Now we have computed the high frequency components of Cnat and
    % Temperature, can perform the PCA analysis on it to get back to kappa_r
end

Z_tmp = zscore(squeeze(high_frequency_tmp));
Z_cnat = zscore(squeeze(high_frequency_Cnat));

X = horzcat(high_frequency_Cnat,high_frequency_tmp);
Z = horzcat(Z_cnat,Z_tmp);
coeff= pca(X);
[~,~,~,~,explainedZ,~] = pca(Z);

if isempty(coeff) | isequal(coeff,eye(2))
    coeff = NaN;
end
% Just send the problematic cases to NaN and then deal with them in the
% following isnan case

if ~isnan (coeff) 
    basis = zeros(2);
    grad = coeff(2,1)./coeff(1,1);
    basis(1,1) = explainedZ(1)./100;
    basis(2,2) = explainedZ(2)./100;
    if n_vargin
		output = grad;
		varargout{1} = sqrt(1 - (basis(4)./basis(1))^2);
       	if verbosity
       	    disp('Output argument 1: Kappa_r');
       	    disp('Output argument 2: Eccentricity');
       	end
else 
    output = NaN;
	for i = 1:n_vargout
    	varargout{i} = NaN;
	end
end

end


function [output] = compute_high_frequencies(input,periodicity,verbosity)

n_cycles = floor(length(input)/periodicity);
if floor(length(input))/periodicity ~= length(input)/periodicity
    n = length(input)/periodicity - floor(length(input))/periodicity;
    if verbosity
        warning(['Input not given for integer number of cycles. Ignoring last ',num2str(n),' data points']);
    end
end

output = NaN(periodicity*n_cycles,1);
		for i = 1:n_cycles
		    output(periodicity*i-periodicity+1:periodicity*i) = input(periodicity*i-periodicity+1:periodicity*i) ...
		        - repmat(nanmean(input(periodicity*i-periodicity+1:periodicity*i)),[periodicity 1]);
		end

end
