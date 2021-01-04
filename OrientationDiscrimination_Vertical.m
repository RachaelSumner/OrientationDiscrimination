% OrientationDescrimination_Vertical

% CREATED:
% Rachael Sumner, December 2020

% EDITED:


% NOTES:

% Runs on Psychtoolbox-3 http://psychtoolbox.org/

% YOU WILL NEED TO ADD YOUR OWN LUMINANCE CORRECTION TO ESTIMATE VOLTAGE NEEDED TO GIVE 44.5 cd/m2
% THE METHOD I USED INCLUDED SCRIPTS THAT ARE NOT MINE TO SHARE.
% THUS VARIABLES mean_lum AND GratingLum ARE UNUSED. SEE GratingGreyLevel.


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     ESSENTIAL PERSONALISATION   %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SavePath = ;%The path to where you want to save data 

ScreenWidth = 29.8;% enter in degrees. Visual angle of screen (degrees). 29.8 assumes sitting 1m away from 53.3cm wide screen
% will be used to ensure stimuli subtend a given degrees of visual angle. 

DeviceID = []; % leave empty for default, else define 

% Personalisation Complete % 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%%%%% BASIC SOUND SETUP %%%%%

InitializePsychSound;
pahandle = PsychPortAudio('Open', DeviceID,[],2,48000,2); 

RLow = MakeBeep(800,.2,48000);
RLow = [RLow; RLow];

RHigh = MakeBeep(1200,.2,48000);
RHigh = [RHigh; RHigh];

WLow = MakeBeep(300,.2,48000);
WLow = [WLow; WLow];

WHigh = MakeBeep(600,.2,48000);
WHigh = [WHigh; WHigh];


%%%%%BASIC SCREEN SETUP

Priority(1);
ListenChar(2); %prevent keyboard input going to the MATLAB window

KbName('UnifyKeyNames');

screens = Screen('Screens'); %For projecting to external screen. Get the number of screens in the computer setup
screenNumber = max(screens); %Can change to 0 for single screen. Otherwise displays on the most external screen (grab max number)

[window, windowRect] = Screen('OpenWindow', screenNumber);
HideCursor;

ScreenRect = Screen ('Rect', window);
pixels_width = ScreenRect(3);
pixels_height = ScreenRect(4);
pixels_per_degree = pixels_width/ScreenWidth;

[xCenter, yCenter] = RectCenter(windowRect); %Finds centre of the screen - Used in Screen('DrawDots',...) for fixation dot

grey = 255 / 2;  

%%%%STIMULI - Luminance corrected

mean_lum = 44.5;

degrees = 4; %degrees of visual angle stimuli will span
num_cycles = 12; %cycles per degree/spatial freq of sine grating 
dia = floor(pixels_per_degree*degrees) - 1; 
radius = dia / 2;

x = 0 :1:dia;

gratingImage_Ref = repmat( ((sin(((x *2* pi * num_cycles)/dia)+pi/2)+1)),[ ( dia + 1 ) 1 ] ); %(x *2* pi * numCycles)/dia freq of sine;  +pi/2 phase; amps scaling
gratingImage_Test1 = repmat( ((cos(((x *2* pi * num_cycles)/dia)+pi/2)+1)),[ ( dia + 1 ) 1 ] ); 
gratingImage_Test2 = repmat( ((-cos(((x *2* pi * num_cycles)/dia)+pi/2)+1)),[ ( dia + 1 ) 1 ] ); 


gratings = whos('gratingImage*');
study_grating= [];

for i = 1: length(gratings)
     
   GratingLum = ((mean_lum).* eval(gratings(i).name));
   
%    GratingGreyLevel = LR.LtoVfun(LR,GratingLum); % This implemented a luminance correction to the grating.

   GratingGreyLevel = gratings(i).name; % in the mean time, this will allow you to see if the rest of the paradigm works but it is not corrected!!!

    % Cut reference grating into a circle
      for y = 0 : dia
            ycentred = y - radius;
            x = sqrt( radius * radius - ycentred * ycentred );
            GratingGreyLevel( y + 1, 1 : floor( radius - x ) ) = grey;
            GratingGreyLevel( y + 1, floor( radius + x + 1 ) : ( dia + 1 ) ) = grey;
      end
          
      study_grating{i} = abs(GratingGreyLevel);

end


%%%% TASK SETUP %%%%

run = 1:44; %size of run 44 = 22 trials per staircase. Must be an even number.
stairs = length(run)/2;

stimulus_on = (0.350);
isi_stairs1 = Shuffle(repmat([0.4 0.5 0.6],1,8));
isi_stairs2 = Shuffle(repmat([0.4 0.5 0.6],1,8));

whichstair = repmat([1 2],1,22); % generate matrix of equal incidence stair 1 and 2
whichstair(3:length(run)) = Shuffle(whichstair(3:length(run))); % keep trial 1 

Reference_Stair1_clockORanti = Shuffle(repmat([-1 1],1,stairs/2));
Rotated_Stair1_clockORanti = Reference_Stair1_clockORanti .* -1;

Reference_Stair2_clockORanti = Shuffle(repmat([-1 1],1,stairs/2));
Rotated_Stair2_clockORanti = Reference_Stair2_clockORanti .* -1;

RandomiseLog = Shuffle(repmat([1:0.5:6],1,2));
GratingDeg = RandomiseLog(1:12);

Adjustment = [0 4 zeros(1,20)+2^0.5];

StartingAngle = 8; 

correctresponsesinarow_1 = 0;
correctresponsesinarow_2 = 0;
j=0;
k=0;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% PARADIGM %%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Participant = input('please enter the study ID: ', 's');
Session = input('please enter session number: ');

Screen ('FillRect', window, grey);

for i = 1:length(run)
    
        Screen('DrawDots', window, [xCenter; yCenter], 10, [0 0 0], [], 2);
        Screen('Flip', window);   
        WaitSecs(0.25)
        
    if whichstair(i) ==1
          j= j + 1;
        %Staircase 1

        % Reference stimulus
      
        if i == 1
            
            x = StartingAngle;
            degangle = x;
            refangle = (degangle/2) * (Reference_Stair1_clockORanti(j));
            
        elseif i>1
            
            if correctresponsesinarow_1 >= 2 || j == 2
                y = log2(x) - log2(Adjustment(j)); %cal octave change
            elseif correctresponsesinarow_1 == 1 && j ~= 2
                y = log2(x);
            elseif correctresponsesinarow_1 == 0 && j ~= 2
                y = log2(x) + log2(Adjustment(j));
            end
            
            by_the_octave = 2^y; % calc log2(x)
            x =  by_the_octave;
            
            degangle =  x; %calc deg = log2(x)
            
            if degangle > 16
                degangle = 16;
            
            elseif degangle < -16
                    degangle = -16;
            end
                
                               
            refangle = (degangle/2) * (Reference_Stair1_clockORanti(j));
                   
        end
        
        
            ReferenceGrating  = Screen('MakeTexture', window, study_grating{1}); 
            Screen('DrawTexture', window, ReferenceGrating, [], [], refangle, [])
            Screen('Flip', window);
            WaitSecs(stimulus_on)    
    
        

        %ISI
        Screen('DrawDots', window, [xCenter; yCenter], 10, [0 0 0], [], 2);
        Screen('Flip', window);   
        WaitSecs(isi_stairs1(j)) 


        % Rotated stim

        rotangle = (degangle/2) * (Rotated_Stair1_clockORanti(j));%rotated angle
        RotatedGrating  = Screen('MakeTexture', window, study_grating{2}); 
        Screen('DrawTexture', window, RotatedGrating, [], [], rotangle)
        
        Screen('Flip', window);


        if Reference_Stair1_clockORanti(j) == -1
            Stair1_results(j).rotation = 'clock'; %
            Stair1_results(j).angle = degangle; %angle size
        elseif Reference_Stair1_clockORanti(j) == 1       
            Stair1_results(j).rotation = 'anti'; %
            Stair1_results(j).angle = degangle; %angle size
        end

        WaitSecs(stimulus_on)      
        
 
        % Respond
        Screen('DrawDots', window, [xCenter; yCenter], 10, [0 0 0], [], 2);
        Screen('Flip', window);
        KbWait;

        [~,~,keyCode]=KbCheck;
        
        if find(keyCode) == KbName('RightArrow')
           
              Stair1_results(j).keypress = 'clock';
   
                    if strfind(Stair1_results(j).keypress, Stair1_results(j).rotation) > 0

                        correctresponsesinarow_1 = correctresponsesinarow_1 +1;
                        
                        PsychPortAudio('FillBuffer', pahandle, RLow);
                        PsychPortAudio('Start', pahandle, 1,0,1);
                        PsychPortAudio('FillBuffer', pahandle, RHigh);
                        PsychPortAudio('Start', pahandle, 1,0,1)

                        Screen('DrawDots', window, [xCenter; yCenter], 10, [0 255 0], [], 2);
                        Screen('Flip', window);
                        WaitSecs(0.2)
                        
                    elseif isempty(strfind(Stair1_results(j).keypress, Stair1_results(j).rotation))

                        correctresponsesinarow_1 = 0;
                        
                        PsychPortAudio('FillBuffer', pahandle, WHigh);
                        PsychPortAudio('Start', pahandle, 1,0,1);
                        PsychPortAudio('FillBuffer', pahandle, WLow);
                        PsychPortAudio('Start', pahandle, 1,0,1)
                    
                        Screen('DrawDots', window, [xCenter; yCenter], 10, [255 0 0], [], 2);
                        Screen('Flip', window);
                        WaitSecs(0.2)
                        
                    end 
                    
                    
                      
          elseif find(keyCode) == KbName ('LeftArrow')  
              
              Stair1_results(j).keypress = 'anti'; 
              
          
                if strfind(Stair1_results(j).keypress, Stair1_results(j).rotation) > 0

                    correctresponsesinarow_1 = correctresponsesinarow_1 +1;

                    PsychPortAudio('FillBuffer', pahandle, RLow);
                    PsychPortAudio('Start', pahandle, 1,0,1);
                    PsychPortAudio('FillBuffer', pahandle, RHigh);
                    PsychPortAudio('Start', pahandle, 1,0,1);

                    Screen('DrawDots', window, [xCenter; yCenter], 10, [0 255 0], [], 2);
                    Screen('Flip', window);
                    WaitSecs(0.2)


                 elseif isempty(strfind(Stair1_results(j).keypress, Stair1_results(j).rotation))

                    correctresponsesinarow_1 = 0;

                    PsychPortAudio('FillBuffer', pahandle, WHigh);
                    PsychPortAudio('Start', pahandle, 1,0,1);
                    PsychPortAudio('FillBuffer', pahandle, WLow);
                    PsychPortAudio('Start', pahandle, 1,0,1)

                    Screen('DrawDots', window, [xCenter; yCenter], 10, [255 0 0], [], 2);
                    Screen('Flip', window);
                    WaitSecs(0.2)

                end 
           
        end   
        



    elseif whichstair(i) ==2 
    
        k = k + 1;
            %Staircase 2

            % Reference stimulus

            if i == 2

                z = StartingAngle;
                degangle = z;
                refangle = (degangle/2) * (Reference_Stair2_clockORanti(k));

            elseif i > 2

                if correctresponsesinarow_2 >= 2 || k == 2
                    y = log2(z) - log2(Adjustment(k)); %cal octave change
                elseif correctresponsesinarow_2 == 1 && k ~= 2
                    y = log2(z);
                elseif correctresponsesinarow_2 == 0 && k ~= 2
                    y = log2(z) + log2(Adjustment(k));
                end

            by_the_octave = 2^y; % calc log2(z)
            z =  by_the_octave;
            
            degangle =  z; %calc deg = log2(z)
                
               if degangle > 16
                    degangle = 16;

                elseif degangle < -16
                        degangle = -16;
                end

                refangle = degangle/2 * (Reference_Stair2_clockORanti(k));

            end
                ReferenceGrating  = Screen('MakeTexture', window, study_grating{1}); 
                Screen('DrawTexture', window, ReferenceGrating, [], [], refangle)
                Screen('Flip', window);
                WaitSecs(stimulus_on)    
        

            %ISI
            Screen('DrawDots', window, [xCenter; yCenter], 10, [0 0 0], [], 2);
            Screen('Flip', window);   
            WaitSecs(isi_stairs2(k)) 


            % Rotated stim

            rotangle = (degangle/2) * (Rotated_Stair2_clockORanti(k));%rotated angle
            RotatedGrating  = Screen('MakeTexture', window, study_grating{3}); 
            Screen('DrawTexture', window, RotatedGrating, [], [], rotangle)
            Screen('Flip', window);


            if Reference_Stair2_clockORanti(k) == -1
                Stair2_results(k).rotation = 'clock'; %
                Stair2_results(k).angle = degangle; %angle size
            elseif Reference_Stair2_clockORanti(k) == 1       
                Stair2_results(k).rotation = 'anti'; %
                Stair2_results(k).angle = degangle; %angle size
            end

            WaitSecs(stimulus_on)


            % Respond
            Screen('DrawDots', window, [xCenter; yCenter], 10, [0 0 0], [], 2);
            Screen('Flip', window);
            KbWait;

            [~,~,keyCode]=KbCheck;
            if find(keyCode) == KbName('RightArrow')
                
                  Stair2_results(k).keypress = 'clock';
                  
                        if strfind(Stair2_results(k).keypress, Stair2_results(k).rotation) > 0

                            correctresponsesinarow_2 = correctresponsesinarow_2 +1;

                            PsychPortAudio('FillBuffer', pahandle, RLow);
                            PsychPortAudio('Start', pahandle, 1,0,1);
                            PsychPortAudio('FillBuffer', pahandle, RHigh);
                            PsychPortAudio('Start', pahandle, 1,0,1)

                            Screen('DrawDots', window, [xCenter; yCenter], 10, [0 255 0], [], 2);
                            Screen('Flip', window);
                            WaitSecs(0.2)

                        elseif isempty(strfind(Stair2_results(k).keypress, Stair2_results(k).rotation))

                            correctresponsesinarow_2 = 0;
                            
                            PsychPortAudio('FillBuffer', pahandle, WHigh);
                            PsychPortAudio('Start', pahandle, 1,0,1);
                            PsychPortAudio('FillBuffer', pahandle, WLow);
                            PsychPortAudio('Start', pahandle, 1,0,1)
                            
                            Screen('DrawDots', window, [xCenter; yCenter], 10, [255 0 0], [], 2)
                            Screen('Flip', window);
                            WaitSecs(0.2)

                        end 

                  
              elseif find(keyCode) == KbName ('LeftArrow')  
                  
                  Stair2_results(k).keypress = 'anti'; 
                                    
                        if strfind(Stair2_results(k).keypress, Stair2_results(k).rotation) > 0

                            correctresponsesinarow_2 = correctresponsesinarow_2 +1;
                                                        

                            PsychPortAudio('FillBuffer', pahandle, RLow);
                            PsychPortAudio('Start', pahandle, 1,0,1);
                            PsychPortAudio('FillBuffer', pahandle, RHigh);
                            PsychPortAudio('Start', pahandle, 1,0,1)

                            Screen('DrawDots', window, [xCenter; yCenter], 10, [0 255 0], [], 2);
                            Screen('Flip', window);
                            WaitSecs(0.2)


                        elseif isempty(strfind(Stair2_results(k).keypress, Stair2_results(k).rotation))

                            correctresponsesinarow_2 = 0;
                            
                            PsychPortAudio('FillBuffer', pahandle, WHigh);
                            PsychPortAudio('Start', pahandle, 1,0,1);
                            PsychPortAudio('FillBuffer', pahandle, WLow);
                            PsychPortAudio('Start', pahandle, 1,0,1)

                            Screen('DrawDots', window, [xCenter; yCenter], 10, [255 0 0], [], 2);
                            Screen('Flip', window);
                            WaitSecs(0.2)

                        end 
                        

            end   


    end
    

    
end

Screen('DrawDots', window, [xCenter; yCenter], 10, [255 0 0], [], 2);
Screen('Flip', window);
WaitSecs(1);

Results.Stair1 = Stair1_results;
Results.Stair2 = Stair2_results;

results_file_name = [SavePath,'OrientationDiscriminationOblique_',num2str(Participant),'_Session',num2str(Session)];
save(results_file_name,'Results')

%%%%END
PsychPortAudio('Close', pahandle)
Screen ('CloseAll');


ShowCursor;
ListenChar (0);
Priority(0)
