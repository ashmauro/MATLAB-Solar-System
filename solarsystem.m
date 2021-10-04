function [p, v] = solarsystem(p, v, mass, stop_time, hide_animation)
% SOLARSYSTEM Solve the gravitational dynamics of n bodies.

if nargin < 5
    hide_animation = false;
end

%% close all open figures
close all;

%% Assign variables to constants

G = 6.673e-11; % Gravitational constant
deltat = 2e3; % Change in time
[n,m] = size(p); % Determine number of variables and store
trail = p; % Variable for body trail

%% Set plot cosmetics

bcolour = {'y', 'b', 'w', 'g', 'r', 'm', 'w'}; % Body Colours
bsize = [100, 35, 15, 30, 45, 30, 35, 20, 10]; % Body sizes

%% Create Plot based on Initial Body data

figure() % Open new figure

for i = 1:n
    
    if m == 2 % 2D Plot
        
        bdata(i) = plot(p(i,1), p(i,2), '.', 'Color', bcolour{i}, 'MarkerSize', bsize(i));% Plot for body
        hold on
        btrail(i) = plot(trail(i,1), trail(i,2), 'w'); % Plot for body trail
        title('2D Solar System Simulation'); % Set 2D title
        
    else % 3D Plot
        
        bdata(i) = plot3(p(i,1), p(i,2), p(i,3),  '.', 'Color', bcolour{i}, 'MarkerSize', bsize(i)); % Plot for body
        hold on
        btrail(i) = plot3(trail(i,1), trail(i,2), trail(i,3), 'w'); % Plot for body trail
        title('3D Solar System Simulation'); %Set 3D title
        
    end
    
end

%Calculate axis limits for variable amount of bodies
axislimit = max(norm(p(1:end)));

if m == 2
    
    % Add extra value on axis limit for a 2 body 2D plot
    axislimit = axislimit + 1e10;
    %2D Labels
    xlabel('Position (m)');
    ylabel('Position (m)');
    %2D Axis limits
    xlim([-axislimit axislimit]);
    ylim([-axislimit axislimit]);
    
else
    %3D Labels
    xlabel('Position (m)');
    ylabel('Position (m)');
    zlabel('Position (m)');
    %3D Axis limits
    xlim([-axislimit axislimit]);
    ylim([-axislimit axislimit]);
    zlim([-axislimit axislimit]);

end

% Set Background Colour
set(gca(),'Color',[0.15 0.15 0.15]);
   

%% Calculating variables for bodies

    for time = 0:deltat:stop_time
        
         for i = 1:n  
             
             F = 0; % Forces are zero
             
                for j = 1:n
                   
                     if i == j

                         continue;

                     end % no self-gravitation
             
                     r = p(j,:) - p(i,:); % Calc of body distances
                     F = F + (G * ((mass(i) * mass(j)) / (norm(r))^3)) * r; % Calc forces acting between bodies
            
                end
            
                a = F / mass(i); % Calc acceleration of each body
                v(i,:) = v(i,:) + (a * deltat); % Calc velocity calculation
                p(i,:) = p(i,:) + (v(i,:) * deltat); % Calc position calculation
                
         end
         
         % Set trail as 3D matrix
         trail(:,:,end+1) = p;

             for q = 1:n

                if m == 2 % 2D Plot
                    
                    %Set 2D plot data
                    set(bdata(q), 'XData', p(q,1), 'YData', p(q,2));
                    set(btrail(q), 'XData', trail(q,1,:), 'YData', trail(q,2,:));

                else % 3D Plot
                    
                    %Set 3D plot data
                    set(bdata(q), 'XData', p(q,1), 'YData', p(q,2), 'ZData', p(q,3));
                    set(btrail(q), 'XData', trail(q,1,:), 'YData', trail(q,2,:), 'ZData', trail(q,3,:));

                end

                drawnow limitrate;

             end    
           
    end
             
end
