clc;
clear;

addpath(genpath('./common'))
% ud_Tri_parameter;

utpr_model_complex_para;

% simple model
% open_system("UTPR_MPC_LPV.slx");

% complex model
open_system("utpr_model_complex.slx")