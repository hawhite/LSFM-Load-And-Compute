%%  LSFM_LoadAndCompute
%   Author: Hamilton White, Postdoctoral Researcher, Brigham and Women's
%   Hospital, Harvard Medical School, and Boston University
%   Copyright (c) 2024 - present: Hamilton White
%
%   License: CC BY-NC-ND
%   (https://creativecommons.org/licenses/by-nc-nd/4.0/)
%
%   Distributed via: https://github.com/hawhite/LSFM-Load-And-Compute
%   Permanent DOI: https://doi.org/10.5281/zenodo.18451522
%
%
%   Core external dependencies:
%       -sigstar
%
%   After the script completes you still need to:
%       -Rename the output file
%       -Have a function to get the variables out of the structures
%       -Compile multiple experiments together in a folder and
%       computationally here in MATLAB
%       -Extract data from the Comps structure, and visualize in a suitable
%       manner
%       -Normalize fonts, size, shape of figures, elements in each panel to
%       the same size/shape/color >> use standards of practice for each
%       journal/review panel as a guide
%
%   General process followed:
%       -Get list of relevant files in the current dir
%       -Per file (saving to Comps structure along the way:
%           >Load in data
%           >Remove flash if existing in the dataset
%           >Normalize data
%           >Suppress some intial datapoints
%           >Differentiate stochastic region of dataset
%           >Compute pairwise correlations of differentiated data
%           >Compute the mean distance moved by each ROI from volume to
%           volume
%           >Compute (if exists) the mean stimulated activation of each
%           neuron
%           >Save outputs to a file
%
%   Temp variables in use:
%       -tmp: raw data input
%       -tmp2: filled, normalized data > input for differentiation
%       -tmp3: temporary output variable > data is then added to the main
%       Comps structure



%% Compute Metrics

% Workspace Cleanup
clc;clear;close all; % Clear command window, close all figure windows, clear variables from workspace

% Start Analysis
disp(sprintf("Starting Analysis: %s",datetime("now"))); % Indicate date and time of analysis start

% Set init variables
skipFirst = 120; % [usual: 120] Timepoints to skip at start of video
dT = 0.5; % [usual: 0.5] Time between volumes (s)
stochasticPeriod = 10; % [usual: 10] Time (min) for the computed period of stochastic activity
flash_period_to_calc = 50; % [usual: 50] Time (timepoints) to calc the mean stimulated activity over
included_differentials = 1; % [usual: 1] Boolean to run the differentiation if desired, else those variables will be set to <empty>
num_parallel_cores = 1; % [usual: 1] Set the number of parallel compute cores to use

% Get file list in current dir
fileList = dir("*.mat"); % Get list of files in current dir
disp(sprintf("Found files:"));
disp(sprintf("\t%s\n",string(extractfield(fileList,"name")'))); % Display that filelist to the user

Comps(size(fileList,1),1) = struct("fn",[],"input",[],"flashPeriods",[],"wholetrace_input",[],"wholetrace_diff",[],"corrs",[],"corrs_tril",[],"roiMovementDistanceMatrix",[],"meanMovementPerNeuronVolume_stochasticRegion",[],"meanStimulatedActivity",[]); % Create empty Comps structure which will contain computational outputs

% Iterate by file in fileList, compute main metrics
parfor (i = 1:size(fileList,1), num_parallel_cores);
    if contains(fileList(i).name,"Experiment")
        continue % Reject anything in the current folder with "Experiment" in the name --> an indicator of analysis output file!
    else
        tmp3 = [];
        disp(sprintf("Now Analyzing: %s, %s",fileList(i).name,datetime("now"))); % Indicate which file is currently being analyzed and the time/date of such analysis being initiated
        tmp = load(fileList(i).name); % Load data
        tmp3.("fn") = string(fileList(i).name); % Capture filename, add to tmp3 structure
        tmp3.("input") = tmp.activity.GCaMP; % Add raw data to tmp3
        [act,flashIDX] = rmFlash(tmp3.input); % Run the rmFalsh function to get a flash-less dataset
        tmp3.("flashPeriods") = flashIDX; % Add the indices of the flash response per video to the tmp3 structure, or empty if no flash exists
        tmp3.("wholetrace_input") = NormalizeFlnrm(act); % Normalize data
        tmp3.wholetrace_input(:,1:skipFirst) = nan; % Suppress 1:skipFirst datapoints due to photobleaching
        %%% Create temporary placeholder variables to maintain size of Comps structure,
        %%% will be replaced below if include_differentials is TRUE
        tmp3.("wholetrace_diff") = [];
        tmp3.("corrs") = [];
        tril_mat = [];
        tmp3.("corrs_tril") = [];
        %%% End temporary placeholder variables
        trace = double(tmp.traces); % Reduce "traces" variable of tracker output from int64 to double
        [K,J] = size(trace,[1,2]); % Get size of trace variable
        Dist_mat = zeros(K,J); % Create empty Distance matrix map
        for k = 1:K
            for j = 2:J
                Dist_mat(k,j) = distCalc(squeeze(trace(k,j,:)),squeeze(trace(k,j-1,:))); % Compute geometric distance between roi and timepoint, enter in corresponding cell of Dist_mat
            end
        end
        tmp3.("roiMovementDistanceMatrix") = Dist_mat; % Adds distance matrix to the tmp3 structure
        tmp3.("meanMovementPerNeuronVolume_stochasticRegion") = mean(Dist_mat(:,(skipFirst+1):(skipFirst+(stochasticPeriod*60./dT))),[1,2]); % Computes mean movement per neuron-volume and adds to tmp3 structure
        if ~isempty(tmp3.flashPeriods)
            tmp3.("meanStimulatedActivity") = meanPeakActivity(tmp3.wholetrace_input,tmp3.flashPeriods,flash_period_to_calc); % Compute mean activation of neuron traces for 50 timepoints after each flash occurs, returns size: [neurons x flashes], add to tmp3 structure
        else
            tmp3.("meanStimulatedActivity") = []; % Else: set empty variable to maintain num of variables across analyses runs
        end
        Comps(i) = tmp3; % Send outputs from temporary structure (tmp3) to main Comps output structure
    end
end

Comps = Comps(~arrayfun(@(x) isempty(x.fn), Comps)); % Remove empty rows of Comps structure before computing differentials, if exist

if included_differentials % Only run if the include_differentials init variable is True, else don't
    WT = struct2mat(3,Comps,[],'wholetrace_input'); % Get whole traces from main structure
    WT_cropped_filled = pagetranspose(fillmissing(pagetranspose(WT(:,(skipFirst+1):(skipFirst+(stochasticPeriod*60./dT)),:)),"nearest",1)); % Crop down data to stochastic region indicated, fill any remaining NaN values in the matrix
    diffMat = pardiffall(WT_cropped_filled, num_parallel_cores); % Run differentiation
    mcorrDiff = mcorr(diffMat); % Compute pairwise correlations by page
    corrs = [];
    for i = 1:size(mcorrDiff,3)
	    tmp=tril(mcorrDiff(:,:,i),-1); tmp(tmp==0)=nan; % Get lower triangle of each page of correlations
	    corrs=safecat(2,corrs,tmp(:)); % concatenate linear array of lower triangle
    end
    parfor (i = 1:size(Comps,1), num_parallel_cores);
        Comps(i).("wholetrace_diff") = diffMat(:,:,i); % Add differentials to Comps structure
        Comps(i).("corrs") = mcorrDiff(:,:,i); % Add correlations to Comps structure
        Comps(i).("corrs_tril") = corrs(:,i); % Add lower triangle linear arrays to Comps structure
    end
end


disp(sprintf("Analysis Complete: %s",datetime("now"))); % Print done


%% Save output
save("Experiment_.mat","Comps","skipFirst","stochasticPeriod","dT","included_differentials"); % Save important output variables

%% Plots of first set (row) of data in Comps set
if included_differentials
    figure(1);clf;
    imagesc(Comps(1).wholetrace_input);

    figure(2);clf;
    edges = linspace(-1,1,50);
    histogram(Comps(1).corrs_tril,'Normalization','probability','BinEdges',edges,DisplayStyle='stairs');
end


%% Functions


%% SafeCat: Safely Concatenate Matrices of Different Sizes

function output = safecat(dim,A,B,blank)
% output = safecat(dim,A,B,blank)
%           like cat function but works with matrices of different size.  
%           Now N-dim arrays supported also.
%           blank is default [NaN]... undefined elements set to blank

    if nargin < 4 blank = NaN; end

    maxdims = max([ndims(A), ndims(B), dim]);
    sizes = [size(A), repmat(1,1,maxdims - ndims(A)); ...
             size(B), repmat(1,1,maxdims - ndims(B))];

    % test for empty matrix
    if isempty(A) 
        output = B;
    elseif isempty(B)
        output = A;
    else
        maxsize = max(sizes);

        if any(sizes(1,:) < maxsize)
            dimexpand = find(sizes(1,:) < maxsize);
            dimexpand = dimexpand(find(dimexpand ~= dim)); % don't expand on 'dim' dimension
            for i = 1:length(dimexpand)
                S.type = '()'; S.subs = {};
                for j = 1:length(maxsize); S.subs{j}=1:maxsize(j); end
                S.subs{dim} = 1:sizes(1,dim);
                S.subs{dimexpand(i)} = (sizes(1,dimexpand(i))+1):maxsize(dimexpand(i));
                A = subsasgn(A,S,blank); 
            end
        end

        if any(sizes(2,:) < maxsize)
            dimexpand = find(sizes(2,:) < maxsize);
            dimexpand = dimexpand(find(dimexpand ~= dim)); % don't expand on 'dim' dimension
            for i = 1:length(dimexpand)
                S.type = '()'; S.subs = {};
                for j = 1:length(maxsize); S.subs{j}=1:maxsize(j); end
                S.subs{dim} = 1:sizes(2,dim);
                S.subs{dimexpand(i)} = (sizes(2,dimexpand(i))+1):maxsize(dimexpand(i));
                B = subsasgn(B,S,blank); 
            end
        end

        output = cat(dim,A,B);
    end
end

%% struct2mat 
%   pulls structure information into matrices
%
%   Example:
%       d = 10; j = 3;
%       for i=1:d rng(2); B(i).data = rand(d,d*2); B(i).idx = linspace(i,i,d)'; end
%       A = struct2mat(1,B,[],{'data'}); % Grab all data
%       IDX = struct2mat(1,B,[1,2,3,5,6],{'idx'}); % Index field for certain index

function output = struct2mat(dim,structure,index,fields)

    if ~iscell(fields) fields = {fields}; end
    if isempty(index) index = 1:numel(structure); end

    output = [];
    for i = 1:length(index)
        substructure = structure(index(i));
        for j = 1:length(fields)
            substructure = getfield(substructure,char(fields(j)));
        end

        output = safecat(dim,output,substructure);
    end
end

%% Compute mean flash period amplitude
function [meanpeak_mat] = meanPeakActivity(data,flash_idx,num_timepoints)
meanpeak_mat = [];
for i = 1:length(flash_idx)
    meanpeak_mat = cat(2,meanpeak_mat,nanmean(data(:,flash_idx(i)+1:flash_idx(i)+(1+num_timepoints)),2));
end
end

%% distCalc: Geometric Distance Calculation for ZYX neuronal centroid data
% Output of neuronal tracker is ZYX, use the following calibration for diSPIM acquisitions:
%       Z = 1 um spacing
%       Y = 0.1625 um spacing
%       X = 0.1625 um spacing
function d = distCalc(firstPt,secondPt)
d = sqrt(1.*(firstPt(1)-secondPt(1))^2 + ...
    0.1625.*(firstPt(2)-secondPt(2))^2 + ...
    0.1625.*(firstPt(3)-secondPt(3))^2);
end

%% rmFlash
function [act,flashIDX] = rmFlash(act)
meanAct = mean(act,1,"omitmissing");
if sum(meanAct,2)==0
    act = act;
    flashIDX = [];
else
    boolIdx = NormalizeFlnrm(meanAct)>50;
    FlashStart = strfind(boolIdx,[0,1])-1;
    flashIDX = [FlashStart',FlashStart'+4];
    for j = 1:size(flashIDX,1)
        act(:,flashIDX(j,1):flashIDX(j,2)) = nan;
    end
end
end

%% TVDiff
function u = tvdiff(data,iter,alph,u0,dx,scale,ep,plotflag,diagflag)
% u = tvdiff(data,iter,alph,dx,u0,scale,ep,plotflag,diagflag);
% Rick Chartrand (rickc@lanl.gov), Apr. 10, 2011
% Please cite Rick Chartrand, "Numerical differentiation of noisy,
% nonsmooth data," ISRN Applied Mathematics, Vol. 2011, Article ID 164564,
% 2011.
%
% Chris Connor (BWH, HMS, BU) - defaults and parameter order made more useful. (1/15/2016)
% Fixed bugs reported at http://math.lanl.gov/~rick/Software/TVDiff/
%
% Inputs:  (First three required; omitting the final N parameters for N < 7
%           or passing in [] results in default values being used.)
%       data        Vector of data to be differentiated.
%
%       iter        Number of iterations to run the main loop.  A stopping
%                   condition based on the norm of the gradient vector g
%                   below would be an easy modification.  No default value.
%
%       alph        Regularization parameter.  This is the main parameter
%                   to fiddle with.  Start by varying by orders of
%                   magnitude until reasonable results are obtained.  A
%                   value to the nearest power of 10 is usally adequate.
%                   No default value.  Higher values increase
%                   regularization strenght and improve conditioning.
%
%       u0          Initialization of the iteration.  Default value is the
%                   naive derivative (without scaling), of appropriate
%                   length (this being different for the two methods).
%                   Although the solution is theoretically independent of
%                   the intialization, a poor choice can exacerbate
%                   conditioning issues when the linear system is solved.
%
%       dx          Grid spacing, used in the definition of the derivative
%                   operators.  Default is the reciprocal of the data size.
%
%       scale       'large' or 'small' (case insensitive).  Default is
%                   'small'.  'small' has somewhat better boundary
%                   behavior, but becomes unwieldly for data larger than
%                   1000 entries or so.  'large' has simpler numerics but
%                   is more efficient for large-scale problems.  'large' is
%                   more readily modified for higher-order derivatives,
%                   since the implicit differentiation matrix is square.
%
%       ep          Parameter for avoiding division by zero.  Default value
%                   is 1e-6.  Results should not be very sensitive to the
%                   value.  Larger values improve conditioning and
%                   therefore speed, while smaller values give more
%                   accurate results with sharper jumps.
%
%       plotflag    Flag whether to display plot at each iteration.
%                   Default is 1 (yes).  Useful, but adds significant
%                   running time.
%
%       diagflag    Flag whether to display diagnostics at each
%                   iteration.  Default is 1 (yes).  Useful for diagnosing
%                   preconditioning problems.  When tolerance is not met,
%                   an early iterate being best is more worrying than a
%                   large relative residual.
%
% Output:
%
%       u           Estimate of the regularized derivative of data.  Due to
%                   different grid assumptions, length( u ) =
%                   length( data ) + 1 if scale = 'small', otherwise
%                   length( u ) = length( data ).

%%% Copyright notice:
% Copyright 2010. Los Alamos National Security, LLC. This material
% was produced under U.S. Government contract DE-AC52-06NA25396 for
% Los Alamos National Laboratory, which is operated by Los Alamos
% National Security, LLC, for the U.S. Department of Energy. The
% Government is granted for, itself and others acting on its
% behalf, a paid-up, nonexclusive, irrevocable worldwide license in
% this material to reproduce, prepare derivative works, and perform
% publicly and display publicly. Beginning five (5) years after
% (March 31, 2011) permission to assert copyright was obtained,
% subject to additional five-year worldwide renewals, the
% Government is granted for itself and others acting on its behalf
% a paid-up, nonexclusive, irrevocable worldwide license in this
% material to reproduce, prepare derivative works, distribute
% copies to the public, perform publicly and display publicly, and
% to permit others to do so. NEITHER THE UNITED STATES NOR THE
% UNITED STATES DEPARTMENT OF ENERGY, NOR LOS ALAMOS NATIONAL
% SECURITY, LLC, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY,
% EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
% RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF
% ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR
% REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED
% RIGHTS.

%%% BSD License notice:
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
%
%      Redistributions of source code must retain the above
%      copyright notice, this list of conditions and the following
%      disclaimer.
%      Redistributions in binary form must reproduce the above
%      copyright notice, this list of conditions and the following
%      disclaimer in the documentation and/or other materials
%      provided with the distribution.
%      Neither the name of Los Alamos National Security nor the names of its
%      contributors may be used to endorse or promote products
%      derived from this software without specific prior written
%      permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
% CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
% INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
% USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
% AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
% LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
% ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

%%% code starts here:
% Make sure we have a column vector.
data = data( : );
% Get the data size.
n = length( data );

% Default checking. (u0 is done separately within each method.)
if nargin < 9 || isempty( diagflag )
    diagflag = 0;
end
if nargin < 8 || isempty( plotflag )
    plotflag = 0;
end
if nargin < 7 || isempty( ep )
    ep = 1e-6;
end
if nargin < 6 || isempty( scale )
    scale = 'small';
end
if nargin < 5 || isempty( dx )
    dx = 1 / n;
end

% Different methods for small- and large-scale problems.
switch lower( scale )

    case 'small'
        % Construct differentiation matrix.
        c = ones( n + 1, 1 ) / dx;
        D = spdiags( [ -c, c ], [ 0, 1 ], n, n + 1 );
        clear c
        DT = D';
        % Construct antidifferentiation operator and its adjoint.
        A = @(x) chop( cumsum( x ) - 0.5 * ( x + x( 1 ) ) ) * dx;
        AT = @(w) ( sum( w ) * ones( n + 1, 1 ) - [ sum( w ) / 2; cumsum( w ) - w / 2 ] ) * dx;
        % Default initialization is naive derivative.
        if nargin < 4 || isempty( u0 )
            u0 = [ 0; diff( data ); 0 ];
        end
        u = u0;
        % Since Au( 0 ) = 0, we need to adjust.
        ofst = data( 1 );
        % Precompute.
        ATb = AT( ofst - data );

        % Main loop.
        for ii = 1 : iter
            % Diagonal matrix of weights, for linearizing E-L equation.
            Q = spdiags( 1 ./ ( sqrt( ( D * u ).^2 + ep ) ), 0, n, n );
            % Linearized diffusion matrix, also approximation of Hessian.
            L = dx * DT * Q * D;
            % Gradient of functional.
            g = AT( A( u ) ) + ATb + alph * L * u;
            % Prepare to solve linear equation.
            tol = 1e-4;
            maxit = 100;
            % Simple preconditioner.
            P = alph * spdiags( spdiags( L, 0 ) + 1, 0, n + 1, n + 1 );
            if diagflag
                s = pcg( @(v) ( alph * L * v + AT( A( v ) ) ), g, tol, maxit, P );
                fprintf( 'iteration %4d: relative change = %.3e, gradient norm = %.3e\n', ii, norm( s ) / norm( u ), norm( g ) );
            else
                [ s, ~ ] = pcg( @(v) ( alph * L * v + AT( A( v ) ) ), g, tol, maxit, P );
            end
            % Update solution.
            u = u - s;
            % Display plot.
            if plotflag
                plot( u, 'ok' ), drawnow;
            end
        end

    case 'large'
        % Construct antidifferentiation operator and its adjoint.
        A = @(v) cumsum(v) * dx;
        AT = @(w) ( sum(w) * ones( length( w ), 1 ) - [ 0; cumsum( w( 1 : end - 1 ) ) ] ) *  dx;
        % Construct differentiation matrix.
        c = ones( n, 1 );
        D = spdiags( [ -c c ], [ 0 1 ], n, n ) / dx;
        D( n, n ) = 0;
        clear c
        DT = D';
        % Since Au( 0 ) = 0, we need to adjust.
        data = data - data( 1 );
        % Default initialization is naive derivative.
        if nargin < 4 || isempty( u0 )
            u0 = [ 0; diff( data ) ];
        end
        u = u0;
        % Precompute.
        ATd = AT( data );

        % Main loop.
        for ii = 1 : iter
            % Diagonal matrix of weights, for linearizing E-L equation.
            Q = spdiags( 1./ sqrt( ( D * u ).^2 +  ep ), 0, n, n );
            % Linearized diffusion matrix, also approximation of Hessian.
            L = DT * Q * D;
            % Gradient of functional.
            g = AT( A( u ) ) - ATd;
            g= g + alph * L * u;
            % Build preconditioner.
            c = cumsum( n : -1 : 1 ).';
            B = alph * L + spdiags( c( end : -1 : 1 ), 0, n, n );
            droptol = 1.0e-2;
            % ---
            % For versions of Matlab > 2013A
            R = ichol( B, struct( 'type', 'ict', 'droptol', droptol ) );
            % Legacy versions use:
            % R = cholinc( B, droptol );
            % ---
            % Prepare to solve linear equation.
            tol = 1.0e-4;
            maxit = 100;
            if diagflag
                s = pcg( @(x) ( alph * L * x + AT( A( x ) ) ), -g, tol, maxit, R', R );
                fprintf( 'iteration %2d: relative change = %.3e, gradient norm = %.3e\n', ii, norm( s ) / norm( u ), norm( g ) );
            else
                [ s, ~ ] = pcg( @(x) ( alph * L * x + AT( A( x ) ) ), -g, tol, maxit, R', R );
            end
            % Update current solution
            u = u + s;
            % Display plot.
            if plotflag
                plot( u, 'ok' ), drawnow;
            end
        end
end
end

% Utility function.
function w = chop( v )
w = v( 2 : end );
end


%% pardiffall Parallelized runtime for tvdiff
function output = pardiffall(act,num_parallel_cores)

% 3D array of neuron traces neuron, frames, epoc(trial) 
% takes the discreate diriviative for each neuron trace 
% outputs the dirivative signals in the same format 

act=double(act);
lthnrm=length(act(:,1,1));
lthfrm=length(act(1,:,1));
lthep=length(act(1,1,:));


parfor (ep=1:lthep, num_parallel_cores);
    
    for nr=1:lthnrm
        tmp=tvdiff(act(nr,:,ep),20,0.1,[],[],'large'); %take deriviative of that trace for that neuron 
        act(nr,:,ep)=tmp(1:lthfrm); 
    end
    ep
end


output=act;    
end

%% MMean: Multi Dimensional Mean Finder
%
%   Examples:
%      mmean(cat(3,[1,1;2,2],[3,3;4,4]),1)
%       1.5000    3.5000
%       1.5000    3.5000
%      mmean(cat(3,[1,1;2,2],[3,3;4,4]),2)
%       1     3
%       2     4
%      mmean(cat(3,[1,1;2,2],[3,3;4,4]),3)
%       2     2
%       3     3

function output = mmean(input,dimlist)

if nargin < 2 dimlist = 1:ndims(input); end

temp = input;
temp(isnan(temp)) = -Inf;

for i = 1:length(dimlist)
    dimswap = 1:ndims(input); dimswap(dimlist(i)) = 1; dimswap(1) = dimlist(i);
    temp = permute(temp,dimswap);
    %temp = sum(temp,dimlist(i));
    temp = mean(temp);
    temp = permute(temp,dimswap);
end
output = squeeze(temp);
end


%% MSum: Multi dimensional Summation

function output = msum(input,dimlist)

if nargin < 2 dimlist = 1:ndims(input); end %#ok<*SEPEX>

temp = input;

for i = 1:length(dimlist)
    dimswap = 1:ndims(input); dimswap(dimlist(i)) = 1; dimswap(1) = dimlist(i);
    temp = permute(temp,dimswap);
    %temp = sum(temp,dimlist(i));
    if size(temp,1) > 1
        temp = nansum(temp);
    end
    temp = permute(temp,dimswap);
end
output = squeeze(temp);
end

%% MStd: Multi Dimensional St.Dev Finder
%
%   Examples:
%      mmean(cat(3,[1,1;2,2],[3,3;4,4]),1)
%       1.5000    3.5000
%       1.5000    3.5000
%      mmean(cat(3,[1,1;2,2],[3,3;4,4]),2)
%       1     3
%       2     4
%      mmean(cat(3,[1,1;2,2],[3,3;4,4]),3)
%       2     2
%       3     3

function output = mstd(input,dimlist)

if nargin < 2 dimlist = 1:ndims(input); end

temp = input;
% temp(isnan(temp)) = -Inf;

for i = 1:length(dimlist)
    dimswap = 1:ndims(input); dimswap(dimlist(i)) = 1; dimswap(1) = dimlist(i);
    temp = permute(temp,dimswap);
    %temp = sum(temp,dimlist(i));
    temp = nanstd(temp);
    temp = permute(temp,dimswap);
end
output = squeeze(temp);
end

%% MCorr: Multi Dimensional Corr Finder
%  Computes corr for each page of a 3D matrix, returning result
%

function [output,tril_vector] = mcorr(input)

sz = size(input,1);
pg = size(input,3);
temp = nan(sz,sz,pg);
tempv = [];

% input(isnan(input)) = -Inf;

input = pagetranspose(input);

for i = 1:pg
    temp(:,:,i) = corr(input(:,:,i),'rows','complete');
    d = temp(:,:,i);
    d(tril(d,-1)==0)=nan;
    tempv(:,i) = d(:);

end
output = temp;
tril_vector = tempv;
end

