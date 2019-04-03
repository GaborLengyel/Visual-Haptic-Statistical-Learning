% Analyses for the study titled: Unimodal statistical learning produces multimodal object-like representations
% The scripts were written by Daniel M. Wolpert, Mate Lengyel, and Gabor Lengyel
% 
% Correspondence to: lengyel.gaabor@gmail.com
%

clear all
clf
close all

load('data.mat') % loading data (can be downloaded from https://osf.io/456qb/)

oldF = cd([pwd,filesep,'myfunctions']);

main = Fig2_fits(Data); % main function that prints the main results showing the average visual and haptic performance, the positive relationship between the two, the comparison across experiments, and the within subject consistency analysis, and generates Figure 2

Fig3_forces(Data) % main function that computes the average pulling forces, the correlations between the pulling and the breakage forces, and generates Figure 3

explicit = Fig4_explicit(Data); % main function that computes the partial correlations after controlling for explicit knowledge and generates Figure 4

% exporting data for the bayes_factors. R script that computes the bayes factors for all of the results
export = nan(max(size(main,1), size(explicit,1)), size(main,2)+size(explicit,2));
export(1:size(main,1),1:size(main,2)) = main;
export(1:size(explicit,1),size(main,2)+1:size(main,2)+size(explicit,2)) = explicit;
cd(oldF)
csvwrite('data_for_R.csv',export)