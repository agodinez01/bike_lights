% Implementation of the Adelson & Bergen (1985) motion energy model.
% Example Matlab code by George Mather, University of Sussex, UK, 2010.
% 
% An expanded version of the code is used by George Mather and Kirsten 
% Challinor as part of a Wellcome Trust funded research project. Initial
% code was partially based on a tutorial forming part of a short course at
% Cold Spring Harbor Laboratories, USA.
% 
% This script is part of an online guide to implementing the Adelson-Bergen
% motion energy model:
%
% http://www.lifesci.sussex.ac.uk/home/George_Mather/EnergyModel.htm
%
% It is free for non-commercial educational use, with appropriate 
% acknowledgement.
%
% The script requires a variable called 'stim' to be loaded in 
% Step 3b. You can use 'AB15.mat' & 'AB16.mat' or input your own stimulus.
%
%--------------------------------------------------------------------------
%           STEP 1: Create component spatiotemporal filters 
%--------------------------------------------------------------------------

clear all; 
close all;

figPath = 'C:\Users\angie\GitRoot\bike_lights\figs\';

scr.measuredFrameRate  = 60; % [Hz] Approximate
params.viewingDistance = 40; % [cm] Approximate
params.framesPerCycle  = scr.measuredFrameRate; % [frames] Equal to the frame rate for 1 Hz

% Step 1a: Define the space axis of the filters
nx = 80;              % Number of spatial samples in the filter
max_x = 2.0;          % Half-width of filter (deg)
dx = (max_x*2)/nx;    % Spatial sampling interval of filter (deg)

% A row vector holding spatial sampling intervals
x_filt = linspace(-max_x,max_x,nx);

% Spatial filter parameters
sx = 0.5;   %standard deviation of Gaussian, in deg.
sf = 1.1;   %spatial frequency of carrier, in cpd

% Spatial filter response
gauss = exp(-x_filt.^2/sx.^2);          % Gaussian envelope
even_x = cos(2*pi*sf*x_filt).*gauss;    % Even Gabor
odd_x = sin(2*pi*sf*x_filt).*gauss;     % Odd Gabor

figure();
plot(even_x, 'b');
hold on;
plot(odd_x, 'g');
xlabel('Space (deg)')
legend('Even Gabor', 'Odd Gabor', box='off')
title('Spatial profile')

saveas(gcf, strcat(figPath, 'Spatial_profile.svg'))

% Step 1b: Define the time axis of the filters
nt = 100;         % Number of temporal samples in the filter
max_t = 0.5;      % Duration of impulse response (sec)
dt = max_t/nt;    % Temporal sampling interval (sec)

% A column vector holding temporal sampling intervals
t_filt = linspace(0, max_t, nt)';

% Temporal filter parameters
k = 100;    % Scales the response into time units
slow_n = 9; % Width of the slow temporal filter
fast_n = 6; % Width of the fast temporal filter
beta =0.9;  % Beta. Represents the weighting of the negative
            % phase of the temporal relative to the positive 
            % phase.

% Temporal filter response (formula as in Adelson & Bergen, 1985, Eq. 1)
slow_t = (k*t_filt).^slow_n .* exp(-k*t_filt).*(1/factorial(slow_n)-beta.*((k*t_filt).^2)/factorial(slow_n+2));
fast_t = (k*t_filt).^fast_n .* exp(-k*t_filt).*(1/factorial(fast_n)-beta.*((k*t_filt).^2)/factorial(fast_n+2));

figure();
plot(slow_t, 'b');
hold on;
plot(fast_t, 'g');
xlabel('Time (sec)')
legend('Slow Filter', 'Fast Filter', box='off')
title('Temporal profile')
saveas(gcf, strcat(figPath, 'Temporal_profile.svg'))

% Step 1c: combine space and time to create spatiotemporal filters
e_slow = slow_t * even_x;   % SE/TS
e_fast = fast_t * even_x;   % SE/TF
o_slow = slow_t * odd_x;    % SO/TS
o_fast = fast_t * odd_x;    % SO/TF

st_filters = {e_slow, e_fast, o_slow, o_fast};

figure()
for i = 1:length(st_filters)
    subplot(1,4,i)
    imagesc(st_filters{i})
    colormap("gray")
    axis off
    axis equal
end

saveas(gcf, strcat(figPath, 'Spatiotemporal_filters.svg'))

%--------------------------------------------------------------------------
%         STEP 2: Create spatiotemporally oriented filters
%--------------------------------------------------------------------------
left_1  = o_fast+e_slow;     % L1
left_2  = -o_slow+e_fast;    % L2
right_1 = -o_fast+e_slow;    % R1
right_2 = o_slow+e_fast;     % R2

st_oriented_filters = {left_1, left_2, right_1, right_2};

figure()
for i = 1:length(st_oriented_filters)
    subplot(1,4,i)
    imagesc(st_oriented_filters{i})
    colormap('gray')
    axis('off')
    axis('equal')
end

saveas(gcf, strcat(figPath, 'Spatiotemporal_oriented_filters.svg'))

%--------------------------------------------------------------------------
%         STEP 3: Convolve the filters with a stimulus
%--------------------------------------------------------------------------

% Step 3a: Define the space and time dimensions of the stimulus

% SPACE: x_stim is a row vector to hold sampled x-positions of the space.
stim_width = 4;  % half width in degrees, gives 8 degrees total
x_stim = -stim_width:dx:round(stim_width-dx);

% TIME: t_stim is a col vector to hold sampled time intervals of the space
stim_dur = 1.5;    %total duration of the stimulus in seconds
t_stim = (0:dt:round(stim_dur-dt))';

% Step 3b Load a stimulus (load 'AB15.mat' or 'AB16.mat') or Create a bitmap
%load('AB15.mat');

bgVal   = 0; % Background color
stimVal = 1; % Stim color

imgSizeX  = 4000;
xs        = 0:imgSizeX;
startXPos = 20;

stimSpeed = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]; % [px/frame]

% how many time intervals to run simulation?
ts = 200;

output = [];
for s=1:length(stimSpeed)
    % figure(s)

    stimBitmap = [];
    for t = 1:ts
        % generat 1D stimulus
        currentT = repmat(bgVal,1,size(xs,2));
        currentX = startXPos + (t-1)*stimSpeed(s);

        if currentX > imgSizeX
            continue
        else
            currentT(currentX) = stimVal;

            stimBitmap(t,:) = currentT;
            % imagesc(stimBitmap);
            % colormap(gray);
            % xlabel('x position');
            % ylabel('time');
        end

    end

    output(s).speed = stimSpeed(s);
    output(s).stim = stimBitmap;

end

% % Make sine wave for flicker
% samples       = 1:ts;
% amplitude     = 0.5;
% phaseShift    = 0;
% verticalShift = 0;
% fullWaveForm = amplitude * sind(360/scr.measuredFrameRate*samples + phaseShift) + verticalShift;
% 
% % Pulse the signal
% for s = 1:length(output)
%     for t = 1:ts
%         if fullWaveForm(t) < 0
%             output(s).stim(t,:) = 0;
%         end
%     end
% end

energy = [];
for s = 1:length(output)
    stim = output(s).stim;

    % Step 3c: convolve
    
    % Rightward responses
    resp_right_1 = conv2(stim, right_1, 'valid');
    resp_right_2 = conv2(stim, right_2, 'valid');
    
    % Leftward responses
    resp_left_1 = conv2(stim, left_1, 'valid');
    resp_left_2 = conv2(stim, left_2, 'valid');
    
    %--------------------------------------------------------------------------
    %         STEP 4: Square the filter output
    %--------------------------------------------------------------------------
    resp_left_1  = resp_left_1.^2;
    resp_left_2  = resp_left_2.^2;
    resp_right_1 = resp_right_1.^2;
    resp_right_2 = resp_right_2.^2;
    
    %--------------------------------------------------------------------------
    %         STEP 5: Normalise the filter output
    %--------------------------------------------------------------------------
    % Calc left and right energy
    energy_right = resp_right_1 + resp_right_2;
    energy_left  = resp_left_1 + resp_left_2;
    
    % Calc total energy
    output(s).total_energy = sum(sum(energy_right))+sum(sum(energy_left));
    
    % Normalise each directional o/p by total output
    RR1 = sum(sum(resp_right_1)) / output(s).total_energy;
    RR2 = sum(sum(resp_right_2)) / output(s).total_energy;
    LR1 = sum(sum(resp_left_1)) / output(s).total_energy;
    LR2 = sum(sum(resp_left_2)) / output(s).total_energy;
    
    %--------------------------------------------------------------------------
    %         STEP 6: Sum the paired filters in each direction
    %--------------------------------------------------------------------------
    output(s).right_Total = RR1 + RR2;
    output(s).left_Total  = LR1 + LR2;
    
    %--------------------------------------------------------------------------
    %         STEP 7: Calculate net energy as the R-L difference
    %--------------------------------------------------------------------------
    output(s).motion_energy = output(s).right_Total - output(s).left_Total;
    energy(:,s) = output(s).motion_energy;
    
    %--------------------------------------------------------------------------
    %         SUPPLEMENTARY CODE: Display summary output and graphics
    %--------------------------------------------------------------------------
    % Display motion energy statistic
    fprintf('\n\nNet motion energy = %g\n\n',output(s).motion_energy);
    
    % Plot the stimulus
    figure(20)
    subplot(3, 4, s);
    imagesc(stim); 
    colormap(gray);
    %axis off
    xlabel('x position');
    ylabel('time');
    title(strcat(num2str(stimSpeed(s)), ' px/frame'));
    %caxis([0, 1.0]);
    %axis equal
    
    % Plot the output:
    %   Generate motion contrast matrix
    energy_opponent = energy_right - energy_left; % L-R difference matrix
    [xv, yv]        = size(energy_left); % Get the size of the response matrix
    energy_flicker  = output(s).total_energy/(xv * yv); % A value for average total energy
    
    % Re-scale (normalize) each pixel in the L-R matrix using average energy.
    motion_contrast = energy_opponent/energy_flicker;
    
    % Plot, scaling by max L or R value
    mc_max = max(max(motion_contrast));
    mc_min = min(min(motion_contrast));
    
    if (abs(mc_max) > abs(mc_min))
        peak = abs(mc_max);
    else
        peak = abs(mc_min);
    end
    
    figure (21)
    subplot(3, 4, s)
    imagesc(motion_contrast); 
    colormap(gray);
    % axis off
    clim([-peak, peak]);
    % axis equal
    title(strcat(num2str(stimSpeed(s)), ' px/frame'));
    %--------------------------------------------------------------------------

end

% % Show quick movie
% figure(); hold on;
% for t =1:ts
%     imagesc(output(3).stim(t,:))
%     colormap(gray);
%     xlabel('x position');
%     ylabel('time')
%     drawnow;
% end


% figure()
% for s = 1:length(output)
%     plot(output(s).right_Total, output(s).left_Total, 'o')
%     hold on;
% end

figure()
for s = 1:length(output)
    plot(output(s).speed, output(s).motion_energy, '.k')
    hold on;
end

f = fit(stimSpeed', energy', 'poly2'); % Fit the data
plot(f, stimSpeed, energy) % Plot the fit

xlabel('stimulus speed (px/frame)');
ylabel('Net motion energy response');