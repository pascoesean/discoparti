%% clear
clear;
clc;

%% load the data
path = '/Users/duyle/Library/CloudStorage/OneDrive-SharedLibraries-MassachusettsInstituteofTechnology/Sean Pascoe - 20.440_final/parti_inputs/pbmcs/';
data = readmatrix([path, 'T Cells/matrix.csv']);
geneNames = readtable([path, 'T Cells/gene_names.csv']);
% clean up dataframe
geneNames(:,1) = [];
geneNames(1,:) = [];
geneNames = table2array(geneNames);

%% make GO matrix
[GOmatrix, GOfullNames, nGenesPerGO, GOcat2Genes] = MakeGOMatrix(data, geneNames, {'MSigDB/c2.cp.v4.0.symbols.gmt', 'MSigDB/c5.all.v4.0.symbols.gmt'}, 10);

%% load in discrete attributes
% In this case, it is just a discrete variable that can values of A, B or C
% for each point.
[discrAttrNames, discrAttr] = ...
    read_enriched_csv([path, 'T Cells/cell_labels.csv'], ',');
% where discrAttr is a matrix of 100 points x 1 attributes. The names of
% the attributes are stored in discrAttrNames.
discrAttrNames = ["T1DM", "Healthy"];

%% parameters
algorithm = 1; % SISAL algorithm
dims = 10; % # of dimensions to analyze
discfeats = ["disease status"];
discmatrix = discrAttr(:, 2);
col = 0; % booleanize all
contfeats = GOfullNames;
contmatrix = GOmatrix;
binsize = 50 / size(data, 1); % pick bin sizes to be 50 cells/bin

%% Run ParTI lite

outputPath = 'Output/TCells/lite';
[arc, arcOrig, pc] = ParTI_lite(data, algorithm, dims, discfeats, discmatrix, col, contfeats, contmatrix, GOcat2Genes, binsize, outputPath);
%writematrix(arc, 'Output/TCells/lite_arc.csv')
%writematrix(arcOrig, 'Output/TCells/lite_arcOrig.csv')
%writematrix(pc, 'Output/TCells/lite_pc.csv')

%% plot the data points
plotParTI(pc,arc,discmatrix)

function plotParTI(pc,arc,annotate_by)
% INPUTS: pc and arc output from parti_lite; column that specifies the
% groups that your rows from pc are in (here, whether or not each cell/row
% came from a T1DM patient or a healthy control. function is currently
% hardcoded to only allow for values of 'T1DM' and something else in this
% array, and will label anything that does not say 'T1DM' as 'Healthy'

figure;
% separate t1d and healthy controls
is_T1DM = strcmp(annotate_by, 'T1DM');

plot3(pc(is_T1DM,1),pc(is_T1DM,2),pc(is_T1DM,3),'.', 'Color', '#606c38','DisplayName', 'T1DM');
hold on;
plot3(pc(~is_T1DM,1),pc(~is_T1DM,2),pc(~is_T1DM,3),'.', 'Color', '#dda15e', 'DisplayName', 'Healthy');

% plot archetypes
plot3(arc(:,1), arc(:,2), arc(:,3), '.r', 'DisplayName', 'Archetype', 'MarkerSize', 15)

% Define faces using indices into arc
% thanks chat
F = [1 2 3 4;  % base
     1 2 5 0;  % sides
     2 3 5 0;
     3 4 5 0;
     4 1 5 0];

for i=1:length(arc(:,1))
    % label archetypes
    text(arc(i,1), arc(i,2), arc(i,3), string("Archetype " + i), 'VerticalAlignment','bottom','HorizontalAlignment','right')
    % make boundary
    face = F(i, F(i,:) > 0);  % skip zeros
    patch('Vertices', arc(:,1:3), 'Faces', face, ...
          'FaceColor', rand(1,3), 'FaceAlpha', 0.05, 'HandleVisibility', 'off')
end
axis equal
legend
box on
xlabel('PC1','fontsize',14);
ylabel('PC2','fontsize',14);
zlabel('PC3','fontsize',14);

hold off;

end


%% Run ParTI
% controls and archetype error estimation.
outputPath = 'Output/ISGhiTCells/full';
[arc, arcOrig, pc, errs, pval] = ParTI(data, algorithm, dims, discfeats, discmatrix, col, contfeats, contmatrix, GOcat2Genes, binsize, outputPath);
writematrix(arc, 'Output/ISGhiTCells/full_arc.csv')
writematrix(arcOrig, 'Output/ISGhiTCells/full_arcOrig.csv')
writematrix(pc, 'Output/ISGhiTCells/full_pc.csv')