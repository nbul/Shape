%% Clear all and initial parameters
clc
clear variables
close all

%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);
%Folders with images
tif8_dir =[filedir, '/tifs_8bit'];
%Folder to save summarised information
mkdir(filedir,'/Summary_shape2');
sum_dir = [filedir, '/Summary_shape2'];

cd(tif8_dir);
files_tif = dir('*.tif');

%% collecting data abount borders and corners
cd(currdir);
collectionv2;

%% Selecting good corners
goodcorners;

%% Centering, rotating and reshaping data
rescale;

%% Fitting and plotting information about borders
borderstat;

%% Fitting and plotting information about corners
cornerstat;

cd(currdir);
close all;
clear variables;
clc





