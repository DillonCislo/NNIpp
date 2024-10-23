%% Scalar Interpolant Test ================================================
% This function tests the Sibson's/Farin's C^1 interpolation method and
% derivative generation functionality of the 'NaturalNeighborInterpolant'
% on a scalar function
%
% by Dillon Cislo 04/26/2020
%==========================================================================

% Uncomment to set
[scriptDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(fullfile(scriptDir, 'mex'));
clear scriptDir

%% Generate Analytic Function and Derivatives =============================
% We test our interpolant on Franke's function, a test case for the
% interpolation of scattered data. You can insert your own favorite
% function (and derivatives) here.
clear; close all; clc;

syms x y
F = franke(x, y);
DF = simplify(gradient(F, [x y]));
DFx = DF(1); DFy = DF(2);

% Convert analytic functions to anonymous functions for evaluation
F = matlabFunction(F, 'Vars', {x, y});
DFx = matlabFunction(DFx, 'Vars', {x, y});
DFy = matlabFunction(DFy, 'Vars', {x, y});

% Create a grid of query points on the domain [0 1] X [0 1]
[X, Y] = meshgrid(0:0.01:1, 0:0.01:1);

trueF = F(X, Y);
trueDFx = DFx(X, Y);
trueDFy = DFy(X, Y);

clear DF x y

%% Set Interpolation Options ==============================================
close all; clc;

nniOptions = struct();
nniOptions.ghostMethod = 'circle';
nniOptions.GPn = 500;
nniOptions.GPr = 2;
nniOptions.GPe = 1;
nniOptions.gradType = 'direct';

switch nniOptions.ghostMethod    
    case 'custom'       
        nniOptions.ghostMethod = 1;       
    case 'circle'        
        nniOptions.ghostMethod = 2;
    case 'edge'        
        nniOptions.ghostMethod = 3;      
end

switch nniOptions.gradType  
    case 'direct'      
        nniOptions.gradType = 1;       
    case 'sibson'       
        nniOptions.gradType = 2;      
end

%% Evaluate Interpolant/Derivatitces ======================================
close all; clc;

% Generate scattered sample points on [0 1] X [0 1]
numPoints  = 2000;
dataPoints = rand([numPoints, 2]);
dataValues = F(dataPoints(:,1), dataPoints(:,2));
dataGrad = [DFx(dataPoints(:,1), dataPoints(:,2)), ...
    DFy(dataPoints(:,1), dataPoints(:,2))];

% Evaluate function
[testF, testDFx, testDFy] = sibsonInterpolantWithGrad( dataPoints(:,1), ...
    dataPoints(:,2), dataValues, nniOptions, X(:), Y(:) );

testF = reshape(testF, size(X));
testDFx = reshape(testDFx, size(X));
testDFy = reshape(testDFy, size(X));

% Report Error Values
FErr = abs(trueF - testF) ./ abs(trueF);
DFxErr = abs(trueDFx - testDFx) ./ abs(trueDFx);
DFyErr = abs(trueDFy - testDFy) ./ abs(trueDFy);

% fprintf('F: Max error = %0.5e\n', max(FErr(:)));
% fprintf('F: RMS Error = %0.5e\n', sqrt(mean(FErr(:).^2)));
fprintf('F: Median Error = %0.5e\n\n', median(FErr(:)));

% fprintf('DFx: Max error = %0.5e\n', max(DFxErr(:)));
% fprintf('DFx: RMS Error = %0.5e\n', sqrt(mean(DFxErr(:).^2)));
fprintf('DFx: Median Error = %0.5e\n\n', median(DFxErr(:)));


% fprintf('DFy: Max error = %0.5e\n', max(DFyErr(:)));
% fprintf('DFy: RMS Error = %0.5e\n', sqrt(mean(DFyErr(:).^2)));
fprintf('DFy: Median Error = %0.5e\n', median(DFyErr(:)));


%% View Results ===========================================================
close all;

fig = figure('Color', 'w', 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

axArray = [];

% Plot True Function Values -----------------------------------------------

axArray(1) = subplot(2,3,1);
surf(X, Y, trueF);
hold on
scatter3(dataPoints(:,1), dataPoints(:,2), dataValues, 'filled', 'r');
hold off
view(145,-2)
xlabel('x')
ylabel('y')
camproj('orthographic')
title('true F')
axis square

axArray(2) = subplot(2,3,2);
surf(X, Y, trueDFx);
hold on
scatter3(dataPoints(:,1), dataPoints(:,2), dataGrad(:,1), 'filled', 'r');
hold off
view(145,-2)
xlabel('x')
ylabel('y')
camproj('orthographic')
title('true \partial F / \partial x')
axis square

axArray(3) = subplot(2,3,3);
surf(X, Y, trueDFy);
hold on
scatter3(dataPoints(:,1), dataPoints(:,2), dataGrad(:,2), 'filled', 'r');
hold off
view(145,-2)
xlabel('x')
ylabel('y')
camproj('orthographic')
title('true \partial F / \partial y')
axis square

% Plot Interpolated Function Values ---------------------------------------

axArray(4) = subplot(2,3,4);
surf(X, Y, testF);
hold on
scatter3(dataPoints(:,1), dataPoints(:,2), dataValues, 'filled', 'r');
hold off
view(145,-2)
xlabel('x')
ylabel('y')
camproj('orthographic')
title('test F')
axis square

axArray(5) = subplot(2,3,5);
surf(X, Y, testDFx);
hold on
scatter3(dataPoints(:,1), dataPoints(:,2), dataGrad(:,1), 'filled', 'r');
hold off
view(145,-2)
xlabel('x')
ylabel('y')
camproj('orthographic')
title('test \partial F / \partial x')
axis square

axArray(6) = subplot(2,3,6);
surf(X, Y, testDFy);
hold on
scatter3(dataPoints(:,1), dataPoints(:,2), dataGrad(:,2), 'filled', 'r');
hold off
view(145,-2)
xlabel('x')
ylabel('y')
camproj('orthographic')
title('test \partial F / \partial y')
axis square

LinkF = linkprop(axArray([1 4]), {'CameraUpVector', 'CameraPosition', ...
    'CameraTarget', 'XLim', 'YLim', 'ZLim'});
LinkDFx = linkprop(axArray([2 5]), {'CameraUpVector', 'CameraPosition', ...
    'CameraTarget', 'XLim', 'YLim', 'ZLim'});
LinkDFy = linkprop(axArray([3 6]), {'CameraUpVector', 'CameraPosition', ...
    'CameraTarget', 'XLim', 'YLim', 'ZLim'});

setappdata(fig, 'StoreTheLink', LinkF);
setappdata(fig, 'StoreTheLink', LinkDFx);
setappdata(fig, 'StoreTheLink', LinkDFy);

