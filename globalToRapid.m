function [varargout] = globalToRapid(varargin)
% A function that takes a global dataset and spits out only the RAPID
% section
nArgs = numel(varargin);
if nArgs
    load('/home/ct/MATLAB/NEMO_runs_downloaded_data/grid_data/section_location_76pt.mat','idx');
else
    error('No inputs given');
end


for i = 1:nArgs
    globalData = varargin{i};
    chop = 0;
    if ndims(globalData) ==2
        globalData = repmat(globalData,[1 1 2]);
        chop = 1;
    end
    if ndims(globalData) == 3
        if size(globalData) == [362,292,64]
            varargout{i} = globalData(idx);
        elseif size(globalData) == [292,362,64]
            fprintf('Lat, Lon appear to be mixed up for input %d. Permuting\n',i)
            globalData = permute(globalData,[2 1 3]);
            varargout{i} = globalData(idx);
        elseif size(globalData,1) == 362 & size(globalData,2) == 292 & size(globalData,3) < 64
            fprintf('Global data not full depth for input %d. Filling in extra levels\n',i)
            fillArray = NaN(362,292,64 - size(globalData,3));
            fillArray = cat(3,globalData,fillArray);
            fillSection = fillArray(idx);
            if chop
                varargout{i} = fillSection(:,1:size(globalData,3)/2);
                fprintf('Output %d: removing extra levels\n',i)
            else
                varargout{i} = fillSection(:,1:size(globalData,3));
                fprintf('Output %d: removing extra levels\n',i)
                
            end
        else
            error('Data is too jumbled (or this function is broken)')
        end
    end
    if ndims(globalData) == 4
        nYears = size(globalData,4);
        fprintf('Assuming fourth dimension is time & selecting rapid section for each of %d time points\n',nYears);
        rapidSec = NaN(76,64,nYears);
        for j = 1:nYears
            rapidSec(:,:,j) = globalToRapid(globalData(:,:,:,j));
        end
        varargout{i} = rapidSec;
    end
            
end
