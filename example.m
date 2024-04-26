addpath('functions/');
addpath('header/');
addpath('data/');

% Please execute the "compile_cpp.m" file to compile C++ functions before using
% the package. Users only need to compile the C++ funciton once.

%% Example 1: 
rng(1);
sim_bm = simulation();
show_simulation(sim_bm);
plot_traj(sim_bm);
plot_intensity(sim_bm.intensity, NaN, NaN, sim_bm.sz);
plot_intensity(sim_bm.intensity, NaN, 10, sim_bm.sz, NaN, true); %10th frame, color image

sim_bm.sim_object = true; % Set this to be TRUE if the intensity profile is from a simulation
sam = SAM(sim_bm);
show_sam(sam);
rng('default');

sim_bm.AIUQ_thr = [0.99, 0.6];
sam = SAM(sim_bm);
show_sam(sam);

sim_bm.index_q_AIUQ = 5:50;
sam = SAM(sim_bm);
show_sam(sam);


%% Example 2:
rng(35);

options = struct();
options.sigma_ou = 4;
options.model_name = "OU";
sim_ou = simulation(options);
show_simulation(sim_ou);

plot_traj(sim_ou);

sim_ou.sim_object = true;
sam_ou = SAM(sim_ou);
show_sam(sam_ou);

rng('default');

plot_MSD(sam_ou, sam_ou.msd_truth);
plot_MSD(sam_ou, sam_ou.msd_truth, NaN, false);

%% Example 3: Simulate BM with smaller frame size and length and get estimated parameters with uncertainty using BM model
rng(1);

options = struct();
options.sz = 100;
options.len_t = 100;
options.sigma_bm = 0.5;
sim_bm = simulation(options);
show_simulation(sim_bm);
rng('default');

plot_traj(sim_bm);

sim_bm.sim_object = true;
sim_bm.uncertainty = true;
sam = SAM(sim_bm);
show_sam(sam);

plot_MSD(sam, sam.msd_truth);
plot_MSD(sam, sam.msd_truth, NaN, false);

%% Example 4
file_path = 'data/intensity_record_BM_beta_0_02_B_40_len_100_frame_size_100.csv';
intensity_ori = readtable(file_path, 'ReadRowNames', true);
intensity = table2array(intensity_ori);

options = struct();
options.mindt = 1;
options.pxsz = 1;
options.sz = [100,100];
options.intensity = intensity;
sam = SAM(options);
show_sam(sam);

options.method = "DDM_fixedAB";
sam = SAM(options);
show_sam(sam);

%% Example 5
tiff_info = imfinfo('data/4pct_PVA_100nm_25C.tif'); 
image_size = size(double(imread('data/4pct_PVA_100nm_25C.tif', 1)));
intensity = zeros(image_size(1), image_size(2), length(tiff_info));

for k = 1:length(tiff_info)
    intensity(:,:,k) = double(imread('data/4pct_PVA_100nm_25C.tif', k));
end

intensity_str = 'SST_array';
plot_intensity(intensity, intensity_str);
plot_intensity(intensity, intensity_str, 15); % 15th frame

options = struct();
options.mindt = 0.0309;
options.pxsz = 0.2930;
options.intensity = intensity;
options.intensity_str = 'SST_array';
options.model_name = "FBM";
options.uncertainty = true;
options.AIUQ_thr = [0.99, 0];
sam = SAM(options);
show_sam(sam);

msd_truth = 0.7589 * sam.d_input;
plot_MSD(sam, msd_truth);


%% Example 6

rng(10);

options = struct();
options.sz = 100;
options.len_t = 100;
options.model_name = "FBM";
sim_fbm = simulation(options);
show_simulation(sim_fbm);
rng('default');

sim_fbm.sim_object = true;
sam_fbm = SAM(sim_fbm);
show_sam(sam_fbm);

sim_fbm.method = "DDM_fixedAB";
sam_fbm = SAM(sim_fbm);
show_sam(sam_fbm);

sim_fbm.method = "DDM_estAB";
sam_fbm = SAM(sim_fbm);
show_sam(sam_fbm);

sim_fbm.method = "DDM_fixedAB";
sim_fbm.index_q_DDM = 3:40;
sam_fbm = SAM(sim_fbm);
show_sam(sam_fbm);

sim_fbm.method = "DDM_estAB";
sim_fbm.index_q_DDM = 3:40;
sam_fbm = SAM(sim_fbm);
show_sam(sam_fbm);


%% Example 7
rng(1);
% Simulation
options = struct();
options.sz = 100;
options.len_t = 100;
options.sigma_bm = 0.5;
sim_bm = simulation(options);
show_simulation(sim_bm);
rng('default');

theta = [2,1];
d_input = 0:10;
model_name = "user_defined";
msd_fn = @user_msd_fn;
MSD_list = get_MSD_with_grad(theta, d_input, model_name, msd_fn);
MSD_list.MSD

% AIUQ method: fitting using user_defined model 
sim_bm.sim_object = true;
sim_bm.user_model_name = "user_defined";
%sim_bm.method = "AIUQ";
%sim_fbm.index_q_DDM = NaN;
sim_bm.num_param = 2;        % For user-defined model, the number of parameters must be given
sim_bm.msd_fn = @user_msd_fn;
sam = SAM(sim_bm);
show_sam(sam);

