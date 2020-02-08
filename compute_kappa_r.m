function [output,varargout] = compute_kappa_r(Cnat_monthly,tmp_monthly,verbosity,varargin)
n_vargin = numel(varargin);

if n_vargin & n_vargin ~= 1
    warning('If entering 4 inputs, must specify scale factor');
    error('Only 3 or 4 outputs can be entered');
elseif n_vargin == 0 & verbosity
    disp('No scale factor specified. Defaulting to determinant');
end

if n_vargin == 1
    if ~strcmpi(varargin{1},'det') & ~strcmpi(varargin{1}, 'ecc') & ~strcmpi(varargin{1}, 'both')
        disp('Scale factor must be eccentricity, determinant or both');
        disp("Enter 'det' for determinant, 'ecc' for eccentricity, 'both' for both")
        error('Scale factor specified incorrectly');
    end
end

Cnat_monthly = squeeze(Cnat_monthly); % Make sure we don't have extraneous dimensions
tmp_monthly = squeeze(tmp_monthly); % Make sure we don't have extraneous dimensions

if ~isvector(tmp_monthly) | ~isvector(Cnat_monthly)
    error('Carbon, Temperature data must be vectors. Error in monthly data');
end

if size(Cnat_monthly,2) ~= 1
    Cnat_monthly = Cnat_monthly'; % Make sure input data is a column vector
end

if size(tmp_monthly,2) ~= 1
    tmp_monthly = tmp_monthly'; % Make sure input data is a column vector
end


if size(tmp_monthly) ~= size(Cnat_monthly)
    error('Carbon, Temperature vectors must be the same size. Error in monthly data.')
end

if verbosity &  isempty(find(~isnan(tmp_monthly))) & isempty(find(~isnan(Cnat_monthly)))
    disp('Temperature and Cnat vectors are all NaN. Assigning NaN output');
    output = NaN;
    varargout{1} = NaN;
    return
end

if nanmean(Cnat_monthly < 50) & nanmean(tmp_monthly > 1500)
   [Cnat_monthly,tmp_monthly] = deal(tmp_monthly,Cnat_monthly);
   if verbosity
       disp("Tmp, Cnat seem to be misplaced. Switching")
   end
end

nan_idx = find(isnan(tmp_monthly) | isnan(Cnat_monthly));
if nan_idx
    tmp_monthly(nan_idx) = 0;
    Cnat_monthly(nan_idx) = 0;
end


high_frequency_tmp = compute_high_frequencies(tmp_monthly,verbosity);
high_frequency_Cnat = compute_high_frequencies(Cnat_monthly,verbosity);


low_frequency_Cnat = Cnat_monthly - high_frequency_Cnat;
low_frequency_tmp = tmp_monthly - high_frequency_tmp;


% In some cases (ie. drift correction), we expect to see a period (ie.
% 30yrs) of temperature and Cnat being zero. Removing these has no effect
% as I've used PCA, but chuck out an error message.

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
        if strcmpi(varargin{1},'det')
            output = grad;
            varargout{1} = 1 - 4*det(basis);
        elseif strcmpi(varargin{1},'ecc')
            output = grad;
            varargout{1} = sqrt(1 - (basis(4)./basis(1))^2);
        elseif strcmpi(varargin{1},'both')
            varargout{1} = 1 - 4*det(basis);
            varargout{2} = sqrt(1 - (basis(4)./basis(1))^2);
            if verbosity
                disp('Output argument 1: Kappa_r');
                disp('Output argument 2: Determinant');
                disp('Output argument 3: Eccentricity');
            end
        end
    else
        output = grad.*(1 - 4*det(basis));
    end
else 
    output = NaN;
    varargout{1} = NaN;
end

end


function [output] = compute_high_frequencies(input,verbosity)

n_years = floor(length(input)/12);
if floor(length(input))/12 ~= length(input)/12
    n = length(input)/12 - floor(length(input))/12;
    if verbosity
        warning(['Input not given for integer number of years. Ignoring last ',num2str(n),' data points']);
    end
end

output = NaN(12*n_years,1);
for i = 1:n_years
    output(12*i-11:12*i) = input(12*i-11:12*i) - repmat(nanmean(input(12*i-11:12*i)),[12 1]);
end

end