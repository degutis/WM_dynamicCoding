%% Define
addpath(genpath('functions'))
addpath(genpath('tdt_3.999C_beta'))
addpath(genpath('spm12'))

global path2data
path2data = '../derivates';
mkdir(path2data)

error('Run manually')

%% Correlation analysis for voxel selection 

% 1: GLMCorr.m
% indicate the experiment: E1 or E2
% to shuffle orientation labels or no: 0/1
% number of iterations
% whether to run feature space smoothing: 0/1
% FWHM of smoothing

GLMCorr('E1',0,1,1,30)
GLMCorr('E1',1,1000,1,30)

GLMCorr('E2',0,1,1,30)
GLMCorr('E2',1,1000,1,30)


% 2: correlationAnalysis.m
% indicate the experiment: E1
% FWHM of smoothing
% whether to conglomerate all IPS into one ROI: 0/1

correlationAnalysis('E1',30,1)
correlationAnalysis('E2',30,1)

%% temporal pSVR

for subject=1:6
    pSVR_temporal('E1',subject,1,1)
    pSVR_gen('E1',subject,1,1)
end

for subject=1:7
    pSVR_temporal('E2',subject,1,1)
    pSVR_gen('E2',subject,1,1)
end


% plot mean temporal 
plotPSVR_temporal_subtraction('E1',1)
plotPSVR_temporal_subtraction('E2',1)

% plot gen 
plotPSVR_gen('E1',1)
plotPSVR_gen('E2',1)


%% PCA subspace analysis
subspace_grandPCA_projection_validated('E1',1,1,1)
subspace_grandPCA_projection_validated('E1',0,8,1)

subspace_grandPCA_projection_validated('E2',1,1,1)
subspace_grandPCA_projection_validated('E2',0,8,1)

%% Simulations
CrossDecodingModel
simulationNoiseAddedRealData("noDist")
simulationNoiseAddedRealData("noiseDist")