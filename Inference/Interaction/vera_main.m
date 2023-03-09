clc 
clear all
close all

TopFolder  = fileparts(fileparts(pwd));
SimFolder  = [TopFolder '\Simulation'];
DataFolder = [SimFolder '\Datasets'];
addpath(([SimFolder, filesep]))
addpath(([fileparts(pwd) '\Diffusion']))
addpath(([fileparts(pwd) '\Allee']))
addpath(([fileparts(pwd) '\Interaction']))


are_we_local = 0;
EXPERIMENT_NR = getenv('SLURM_ARRAY_TASK_ID'); 

if are_we_local == 1
    df = dir(DataFolder);
    df = df(3:end);
    load([DataFolder '\' df(end).name])
elseif are_we_local == 0
    load('20230226_1743_dataset_real.mat')
end

inference_main(are_we_local,EXPERIMENT_NR,observed_cells)