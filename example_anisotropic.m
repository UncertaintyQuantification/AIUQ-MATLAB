addpath('functions/');
addpath('header/');


%% Example 1: Simulate anisotropic BM and get estimated parameters using anisotropic BM model
rng(11);

options = struct();
options.sigma_bm = [0.5, 0.1];
sim_bm = aniso_simulation(options);
show_aniso_simulation(sim_bm);
plot_traj(sim_bm);

plot_intensity(sim_bm.intensity, NaN, NaN, sim_bm.sz);
plot_intensity(sim_bm.intensity, NaN, 10, sim_bm.sz, NaN, true); %10th frame, color image

sim_bm.sim_object = true; % Set this to be TRUE if the intensity profile is from a simulation
%sim_bm.AIUQ_thr = [0.999, 0];
aniso_sam = aniso_SAM(sim_bm);
show_aniso_sam(aniso_sam);
rng('default');
% Plot estimated MSD
plot_MSD(aniso_sam, aniso_sam.msd_truth);


%%  Example 2: Simulate isotropic BM and get estimated parameters using (an)isotropic BM models
rng(25);

options = struct();
options.sigma_bm = 0.5;
options.sz = 100;
options.len_t = 100;
sim_bm = simulation(options);
show_simulation(sim_bm);
plot_traj(sim_bm);
plot_intensity(sim_bm.intensity, NaN, NaN, sim_bm.sz);

sim_bm.sim_object = true; % Set this to be TRUE if the intensity profile is from a simulation
%sim_bm.AIUQ_thr = [0.999, 0];

% AIUQ method: use anisotropic BM as fitted model 
aniso_sam = aniso_SAM(sim_bm);
show_aniso_sam(aniso_sam);
% AIUQ method: use isotropic BM as fitted model
sam = SAM(sim_bm); 
show_sam(sam);
rng('default');
% Plot estimated MSD
plot_MSD(aniso_sam, aniso_sam.msd_truth);
plot_MSD(sam, sam.msd_truth);



%% Example 3: Simulate anisotropic OU and get estimated parameters using anisotropic OU model with uncertainty
rng(35);

options = struct();
options.model_name = "OU";
options.sz = 100;
options.len_t = 100;
options.sigma_ou = [4, 1.5];
options.rho = [0.95, 0.85];
options.sigma_p = 1;

aniso_sim_ou = aniso_simulation(options);
show_simulation(aniso_sim_ou);
% Plot particle trajectory and intensity profile for time 1
plot_traj(aniso_sim_ou);
plot_intensity(aniso_sim_ou.intensity, NaN, NaN, aniso_sim_ou.sz);

aniso_sim_ou.sim_object = true; 
%aniso_sim_ou.AIUQ_thr = [0.999, 0];
%aniso_sim_ou.uncertainty = true;
% AIUQ method: use anisotropic OU as fitted model 
aniso_sam = aniso_SAM(aniso_sim_ou);
show_aniso_sam(aniso_sam);
rng('default');
% Plot estimated MSD
plot_MSD(aniso_sam, aniso_sam.msd_truth);





%% Example 4: Simulate anisotropic FBM and get estimated parameters with uncertainty
rng(36);

options = struct();
options.model_name = "FBM";
options.sigma_fbm = [0.5, 1];
options.H = [0.25, 0.45];

aniso_sim_fbm = aniso_simulation(options);
show_simulation(aniso_sim_fbm);
% Plot particle trajectory and intensity profile for time 1
plot_traj(aniso_sim_fbm);
plot_intensity(aniso_sim_fbm.intensity, NaN, NaN, aniso_sim_fbm.sz);

aniso_sim_fbm.sim_object = true; 
aniso_sim_fbm.AIUQ_thr = [0.999, 0];
aniso_sim_fbm.uncertainty = true;
% AIUQ method: use anisotropic FBM as fitted model 
aniso_fbm_sam = aniso_SAM(aniso_sim_fbm);
show_aniso_sam(aniso_fbm_sam);
rng('default');
% Plot estimated MSD
plot_MSD(aniso_fbm_sam, aniso_fbm_sam.msd_truth);



%% Example 5: Load csv file and estimate parameters using anisotropic FBM model
file_path = 'data/intensity_record_FBM_anisotropic_sigma_0.4_0.8_H_0.25_0.45_len_100_frame_size_100.csv';
intensity_ori = readtable(file_path, 'ReadRowNames', true);
intensity = table2array(intensity_ori);

options = struct();
options.mindt = 1;
options.pxsz = 1;
options.intensity = intensity;
options.uncertainty = true;
options.sz = [100,100];
options.model_name = "FBM";

aniso_sam = aniso_SAM(options);
show_aniso_sam(aniso_sam);

% Plot estimated MSD in x,y directions. First construct true/reference if exist, as a dimension delta t*2 matrix
true_param = [0.64, 0.16; 0.9, 0.5];
msd_truth = zeros(aniso_sam.len_t, size(true_param, 2));
for i = 1:size(true_param, 2)
    msd_truth(:,i) = get_MSD(true_param(:, i), aniso_sam.d_input, aniso_sam.model_name);
end

plot_MSD(aniso_sam, msd_truth);




%% Example 6: Simulate anisotropic OU+FBM and get estimated parameters using anisotropic OU+FBM model
rng(8);

options = struct();
options.model_name = "OU+FBM";
options.sigma_fbm = [2, 1];
options.sigma_ou = [1.5, 1];
options.H = [0.3, 0.225];
options.rho =[0.85, 0.8];

aniso_sim_ou_fbm = aniso_simulation(options);
show_simulation(aniso_sim_ou_fbm);
rng('default');
% Plot particle trajectory and intensity profile for time 1
plot_traj(aniso_sim_ou_fbm);
plot_intensity(aniso_sim_ou_fbm.intensity, NaN, NaN, aniso_sim_ou_fbm.sz);

aniso_sim_ou_fbm.sim_object = true; 
aniso_sim_ou_fbm.AIUQ_thr = [0.999, 0];
aniso_sim_ou_fbm.uncertainty = true;
aniso_ou_fbm_sam = aniso_SAM(aniso_sim_ou_fbm);
show_aniso_sam(aniso_ou_fbm_sam);
% Plot estimated MSD
plot_MSD(aniso_ou_fbm_sam, aniso_ou_fbm_sam.msd_truth);



%% Example 7: User defined MSD structure
rng(11);

options = struct();
options.sigma_bm = [0.5, 0.1];
options.sz = 100;
options.len_t = 100;
sim_bm = aniso_simulation(options);
show_aniso_simulation(sim_bm);

% show MSD and MSD gradient with a simple example 
theta_x = [2, 1];
d_input = 0:10;
model_name = 'user_defined';
MSD_list_x = get_MSD_with_grad(theta_x, d_input, model_name, @user_msd_fn);
disp(MSD_list_x.MSD);

sim_bm.sim_object = true;
sim_bm.user_model_name = "user_defined";
sim_bm.num_param = 2; % For user-defined model, the number of parameters must be given
sim_bm.msd_fn = @user_msd_fn;
% AIUQ method: fitting using user_defined model 
aniso_sam = aniso_SAM(sim_bm);
show_aniso_sam(aniso_sam);
% Plot estimated MSD
plot_MSD(aniso_sam, aniso_sam.msd_truth);




%% Example 8: Load tif file (SST_array) from real experiment and est parameters using FBM model
tiff_info = imfinfo('data/16pct_DSCG_200nm_3.tif'); 
image_size = size(double(imread('data/16pct_DSCG_200nm_3.tif', 1)));
intensity = zeros(image_size(1), image_size(2), length(tiff_info));

for k = 1:length(tiff_info)
    intensity(:,:,k) = double(imread('data/16pct_DSCG_200nm_3.tif', k));
end

intensity_str = 'SST_array';
plot_intensity(intensity, intensity_str);
plot_intensity(intensity, intensity_str, 100); % 100th frame

%%% FBM
options = struct();
options.mindt = 0.0618;
options.pxsz = 0.29;
options.intensity = intensity;
options.intensity_str = 'SST_array';
options.model_name = "FBM";
options.uncertainty = true;
options.AIUQ_thr = [0.999, 0];
% sam_DSCG = aniso_SAM(options);
% show_aniso_sam(sam_DSCG);

% msd_truth_x = 0.0167 * sam_DSCG.d_input;
% msd_truth_y = 0.0027 * sam_DSCG.d_input;
% msd_truth = [msd_truth_x', msd_truth_y'];
% MPT
MSD_MPT_x = readmatrix("data/16pct_DSCG_200nm_3_MSDMicx.csv");
MSD_MPT_y = readmatrix("data/16pct_DSCG_200nm_3_MSDMicy.csv");

% plot_MSD(sam_DSCG, msd_truth);
% hold on;
% plot(log10(sam_DSCG.d_input(2:end)), log10(MSD_MPT_x), '--', 'Color', '#D95319', 'LineWidth', 2);
% plot(log10(sam_DSCG.d_input(2:end)), log10(MSD_MPT_y), '-', 'Color', '#D95319', 'LineWidth', 2);


options.AIUQ_thr = [0.99, 0];
sam_DSCG = aniso_SAM(options);
show_aniso_sam(sam_DSCG);

msd_truth_x = 0.0167 * sam_DSCG.d_input;
msd_truth_y = 0.0027 * sam_DSCG.d_input;
msd_truth = [msd_truth_x', msd_truth_y'];
plot_MSD(sam_DSCG, msd_truth);
hold on;
plot(log10(sam_DSCG.d_input(2:length(MSD_MPT_x)+1)), log10(MSD_MPT_x), '--', 'Color', '#D95319', 'LineWidth', 2, 'HandleVisibility', 'on');
plot(log10(sam_DSCG.d_input(2:length(MSD_MPT_y)+1)), log10(MSD_MPT_y), '-', 'Color', '#D95319', 'LineWidth', 2, 'HandleVisibility', 'on');
legend('Reference x','Reference y','SAM x','SAM y', 'MPT x', 'MPT y', 'Location', 'NorthWest', 'FontSize', 11);
 