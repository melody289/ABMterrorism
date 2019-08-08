% This is the Set up for the Model

%feat is the population demographics for the country of interest.
% The features are in the following order
% 1) Percent of 15 and older which are married, I found this from UN
% 2) Shape of gamma of wealth, I have choosen 2 for france and 1 for Libya
% 3) Percent of population which is NOT islamic
% 4) illiterate percent
% 5) No education, literate
% 6) some primary education
% 7) some secondary education
% 8) some tertiary, Note that 4-8 should sum to 1
%9) exposure to crime intialization
% 10) Number of Perp


% This is for France
%feat = [ 0.45, 2, 0.95,0, 0.015, 0.164, 0.588, 0.233, 0, 2];
% This is for Libya
feat = [ 0.43, 1, 0.005,0.086, 0.221, 0.252, 0.2932, 0.2338, 1, 0];
% This is for Pakistan
%feat = [ 0.63, 0.5, 0.04,0.086, 0.173, 0.172, 0.317, 0.057, 1,0];

% This outputs the Network A2 (F in graph form) and characteristic matrix
% cha, total population n
[A2, cha, F, n] = popDem(feat);


% This is the power level need for a terrorist to be able to recruit or
% something
powL = 15;
powL2 = 7;
powLf = 10;
arrest = 20;

% This is the scale for the logistic distribution, I choose it to be sharp
sL = 0.1;
sF = 0.1;
sP = 0.1;

% I want this to be a gentlier rise
sT = 100;

% This is the means for the logistic distribution values choosen so mean is
% zero, but can be changes
muL = 0; muP =-1; muF =0;

cha2 = cha;
% This makes a successful run if wT*[finance power) >56
muT = 0;



% This is the weights for equal spread
% Code used, but need to use same values each time and I did not record
% seed, so I will upload it
%ww = lhsdesign(samp, 8*3 + 3);
%ww(:, 1:24) = ww(:, 1:24).*2 -1;
load('weightsforall.mat')

% This is the ones which are weighted towards our hypothesis

load('weightsspecific.mat')

% putting weights together
ww = [ ww;ww2];

ttday = zeros(50,7);
big = zeros(50,7,40);

%%

% jjj is the specific weight choosen

for jjj = [15, 31, 35, 40]
    
    tic
wL = ww(jjj,1).*ww(jjj,2:8);
wF = ww(jjj,9).*ww(jjj,10:16);
wP = ww(jjj,17).*ww(jjj,18:23);
wT = ww(jjj,25:27);
parfor k = 1:50

    A = A2;
  cha = cha2;
% % % Weights for Leader
% % wL = 0.5.*[ 0.5    0.5   0  0.1    0.5    0.5    0.5];
% % 
% % % Weights for Financier
% % 
% % wF = 0.5.*[-0.5 0.5 2 0 0.5  0 1  ];
% % 
% % % Weights for perpetrator
% % 
% % wP = 2.*[ 0  -0.5 -0.5 0.5  1  -0.5];
% % 
% % % Weight for terror attack
% % 
% % wT = [ 0.5,0.5, 0.5];

% counts for suc, even, civial kill, terrorist kill

suc =0;
even =0;
nkill = 0;
tkill =0;

wta = 3;
% Number of connections (neighbor, friend)




%
%%%%%%%%%%%%%%%%
%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%%%%%%%%
%%%%%%%%%%%%%%%
%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%%%%%%%%
it = 4000;
%

for jj = 1:it

% One round, we randomly select ~ 10%

s = RandStream('mt19937ar','Seed',2510, 'NormalTransform', 'Ziggurat');
RandStream.setGlobalStream(s)

select = randperm(n, round(.1*n) +1);

zl = length(select);

for i = 1:zl
    % zi is the initiator and zs is the one interacted with
    zi = select(i);
    % This finds one of the neighbors of select(i)
    zs = find(A(zi,:));
    
    
   
    if(~isempty(zs))
    zs = zs(randi([1, length(zs)],1));
     u = rand;
     
     % What happens when one group type to another interact

%%%
%%%
%%%

     % A civilian is choosen
    if(cha(zi,1) == 0)
        
        % with Military
        if(cha(zs,1) == 1)
            % Increase or decrease bias towards military
        cha( zi , 11) = cha(zi, 11) + (u/5)*(-1)^(cha(zi, 11) < u);
         
         % with terrorist
        elseif( cha(zs,1) ~= 0)  
            
            % work
            if(powL*3 < cha(zs,9) && cha(zi,10) < .2 )
                 if(1/(1+exp(-(wP*[cha(zs,2:7)']+.5 - muP)/sP)) > u && cha(zs, 13) < 1)
                        
                        cha(zs,10) = 0; 
                        cha(zi,9) = cha(zi,9) + 1;
                        
                        % connecting to cell leader
                        z2 = intersect(find(A(zi,:)), find(cha(:,1) == 4, 1));
                        A(z2,zs) = 200;
                        A(zs, z2) = 10;
                        
                        % connecting to cell
                        z2 = intersect(find(A(zi,:)), find(cha(:,1) == 2));
                        z2 =  z2(A(z2, zi)==50);
                        A(z2,zs) = 50;
                        A(zs, z2) = 50;
                        
                       
                        cha(zs,1) = 2;
                 end
            else
%             % Increase/decrease bias towards terrorist
        cha(zi, 10) = cha(zi, 10) + cha(zs,1)*(u/15)*(-1)^(cha(zi, 10) < u);
            end
        
        end
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%
%%%       
        
      % Foot Soldier  
    elseif(cha(zi,1) == 2)
      if(cha(zi,13) >0)
        cha(zi,13)= 0;
%%%
%%%
%%%        
      else
        %FT to C
         if(cha(zs,1) == 0 )
             
             % If there is enough power and low bias, try to recuit
           if(cha(zi, 9) > powL2 && cha(zs,10) < .2)
                % This is the logistic distribution for each type of
                % Terrorist (L, F, P)
                
                a2 = [1/(1+exp(-(wP*cha(zs,2:7)' - muP)/sP)),  1/(1+exp(-(wF*[cha(zs,2:7)'; (cha(zs,8) ==2)] - muF)/sF)), 1/(1+exp(-(wL*[cha(zs,2:7)'; (cha(zs,8) ==1)] - muL)/sL))];
                 
                % Becoming  FS
                if(a2(1)> 0.5)
                        
                        cha(zs,10) = 0; 
                        cha(zi,9) = cha(zi,9) + 1;
                        
                          % connecting to cell
                        
                        z2 = intersect(find(A(zi,:)== 50), find(cha(:,1) == 2));
                        
                     
                        A(z2,zs) = 50;
                        A(zs, z2) = 50;
                        
                        
                       
                        % connecting to cell leader
                        z2 = intersect(find(A(zi,:)== 50), find(cha(:,1) == 4));
                        %z2 =  find(A(z2, zi)==200, 1);
                        A(z2,zs) = 200;
                        A(zs, z2) = 50;
                        
                      
                        cha(zs,1) = 2;
                        
                        
                        
                elseif(a2(2) > 0.5 || a2(3) > 0.5)
                    
                    
                    % Becoming financier
                    if(a2(2) >  a2(3) && cha(zs,4) > 0 )
                        cha(zs,1) = 3;
                        cha(zs,10) = 0; 
                        cha(zi,9) = cha(zi,9) + 1;
                        
                         % connecting to cell leader
                        z2 = intersect(find(A(zi,:)==50), find(cha(:,1) == 4));
                        A(z2,zs) = 201;
                        A(zs, z2) = 201;
                        
                       
                        % Becoming cell leader
                    else
                        z2 = intersect(find(A(zi,:)==50), find(cha(:,1) == 4));
                        z3 = intersect(find(A(z2,:)== 50), find(cha(:,1) == 2));
                        zp = length(z3);
                        if(zp > 5)
                            cha(zs,1) = 4;
                            cha(zs,10) = 0; 
                            cha(zi,9) = cha(zi,9) + 1;
                            % connect to other cell leader
                            A(z2,zs) = 201;
                            A(zs, z2) = 201;
                            
                            % Dividing cell
                            z4 = randperm(zp, 3);
                            z4 = z3(z4);
                            A(z2,z4) = 0;
                            A(z4, z2) = 0;
                            A(z4,zs) = 200;
                            A(zs, z4) = 50;
                            
                            z5 = setdiff(z3, z4);
                            A(z4, z5 ) = 0;
                            A(z5, z4 ) = 0;
                            A(z4, z4) = 50;
                            A(z5, z5) = 50;
                            
                           A =  A - diag(diag(A));
                           
                        else
                            cha(zs, 10) = cha(zs, 10) + rand/5;  
                        end
                        
                    end
                    
                else
                cha(zs,10) =  cha(zs, 10) -rand/5;
                cha(zi,9) = cha(zi,9) -1;
                end
            else
             cha(zs, 10) = cha(zs, 10) + cha(zs,1)*(u/15)*(-1)^(cha(zs, 10) < u);    
           end
           
           %
           % % Meets another terrorist
         elseif( cha(zs,1) ~= 1)
          cha(zi,9) = cha(zi,9) + 1;   
          cha(zs,9) = cha(zs,9) + 1;  
%%%
%%%
%%%
          % m hide
         elseif(cha(zs,9) > cha(zi,9))
             cha(zi,13) =  1;
         end
     
      end  
%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%%%%%
%%%%      
    
         %Financier
    elseif(cha(zi,1) == 3)
  if(cha(zi,13) >0) 
      
       
      
       cha(zi,13) = 0;
%%%%%%%%%
%%%
  else
           
%%%
      
        %meets Civilian
         if(cha(zs,1) == 0 )
             
             % If there is enough power and low bias, try to recuit
            if(cha(zi, 9) > powL2 && cha(zs,10) < .2)
               % This is the logistic distribution for each type of
                % Terrorist (L, F, P)
                a2 = [1/(1+exp(-(wP*[cha(zs,2:7)'] - muP)/sP)),  1/(1+exp(-(wF*[cha(zs,2:7)'; (cha(zs,8) ==2)] - muF)/sF)), 1/(1+exp(-(wL*[cha(zs,2:7)'; (cha(zs,8) ==1)] - muL)/sL))];
                  
                if(a2(1)> u)
                        
                        cha(zs,10) = 0; 
                        cha(zi,9) = cha(zi,9) + 1;
                        
                        % connecting to cell leader
                        z2 = intersect(find(A(zi,:)==201), find(cha(:,1) == 4));
                        z2 = z2(end);
                        A(z2,zs) = 200;
                        A(zs, z2) = 50;
                        
                        % connecting to cell
                        z3 =  find(A(z2, :)==200);
                        A(z3,zs) = 50;
                        A(zs, z3) = 50;
                        
                       
                        cha(zs,1) = 2;
                        
                        
                elseif(a2(2) > u || a2(3) > u)
                    if(a2(2) >  a2(3) && cha(zs,4) > 0)
                        cha(zs,1) = 3;
                        cha(zs,10) = 0; 
                        cha(zi,9) = cha(zi,9) + 1;
                        
                         % connecting to cell leader
                        z2 = intersect(find(A(zi,:)), find(cha(:,1) == 4));
                        z2 = z2(end);
                        A(z2,zs) = 201;
                        A(zs, z2) = 201;
                        
                        
                    else
                        z2 = intersect(find(A(zi,:)), find(cha(:,1) == 4));
                        z3 = zeros(length(z2),1);
                        for i2 = 1:length(z2) 
                            z3(i2) = length( find(A(z2(i2), :)==200));
                        end
                        [zp, b] = max(z3);
                         
                        if(zp > 5)
                            z2 = z2(b);
                            z3 = find(A(z2,: )== 200);
                            cha(zs,1) = 4;
                            cha(zs,10) = 0; 
                            cha(zi,9) = cha(zi,9) + 1;
                            
                            % connect to other cell leader
                            A(z2,zs) = 201;
                            A(zs, z2) = 201;
                           
                            
                             % Dividing cell
                            z4 = randperm(zp, 3);
                            z4 = z3(z4);
                            A(z2,z4) = 0;
                            A(z4, z2) = 0;
                            A(z4,zs) = 200;
                            A(zs, z4) = 50;
                            
                            z5 = setdiff(z3, z4);
                            A(z4, z5 ) = 0;
                            A(z5, z4 ) = 0;
                            A(z4, z4) = 50;
                            A(z5, z5) = 50;
                            
                           A =  A - diag(diag(A));
                            
                        else
                            cha(zs, 10) = cha(zs, 10) + rand/5;  
                        end
                        
                    end
                            
                else
                cha(zs,10) =  cha(zs, 10) -rand/5;
                cha(zi,9) = cha(zi,9) -1;
                end
            else
             cha(zs, 10) = cha(zs, 10) + cha(zs,1)*(u/15)*(-1)^(cha(zs, 10) < u);    
                       
            end
 %%%%
 %%%%
 %%%%
 
            % Financier meet FS
         elseif( cha(zs,1) == 2)
          cha(zi,9) = cha(zi,9) + 1;   
          cha(zs,9) = cha(zs,9) + 1; 
 %%%%
 %%%%
 %%%        
          
          %Financier meets L
         elseif( cha(zs,1) == 4)
             z2 = find(A(zs,:) ==200);
             if(cha(zi,9) >powLf && sum(cha(z2 , 1) ==2) > 2)
                 
                 % %
                 cha(zs, 9) = cha(zs, 9) + 2;
                 
                 
                 % Starts the cycle, finances cell
                 cha([z2,zs],12) = cha([z2,zs],12) +  2;
                        
             else
             cha(zi,9) = cha(zi,9) + 1;   
             cha(zs,9) = cha(zs,9) + 1; 
             end
  %%%%
 %%%%
 %%%       Financier meet Military            
         elseif(cha(zs,1) ==1 && (cha(zs,9)  > cha(zi,9) + 2))   
             cha(zi,13) =   1;
             
         end
  end
 %%%%
 %%%%%
 %%%%
 %%%%%%%%%%%%%%%%%%%%%%%%
 %%%%
   
 
        %Cell leader
    elseif(cha(zi,1) == 4)
     
        if(cha(zi,13) >0)
            cha(zi, 13) = 0;
        else
               
        %meets civilian
         if(cha(zs,1) == 0 )
             
             % If there is enough power and low bias, try to recuit
            if(cha(zi, 9) > powL2 && cha(zs,10) < .2)
                % This is the logistic distribution for each type of
                % Terrorist (L, F, P)
                a2 = [1/(1+exp(-(wP*cha(zs,2:7)' - muP)/sP)),  1/(1+exp(-(wF*[cha(zs,2:7)'; (cha(zs,8) ==2)] - muF)/sF)), 1/(1+exp(-(wL*[cha(zs,2:7)'; (cha(zs,8) ==1)] - muL)/sL))];
                  
                if(a2(1)> u)
                        
                        
                        cha(zs,10) = 0; 
                        cha(zi,9) = cha(zi,9) + 1;
                        
                        % connecting to cell leader with correct weights
                       
                        A(zi,zs) = 200;
                        A(zs, zi) = 50;
                        
                        % connecting to cell
                        z2 = intersect(find(A(zi,:)==200), find(cha(:,1) == 2));
                       
                       
                        A(z2,zs) = 50;
                        A(zs, z2) = 50;
                        cha(zs,1) = 2;
                        if(length(A) ~= n)
                         msg = 'Error occurred!!! Annoying 755';
                          error(msg)
                         end
                        
                elseif(a2(2) > u || a2(3) > u)
                    if(a2(2) >  a2(3)&& cha(zs,4) > 0)
                        cha(zs,1) = 3;
                        cha(zs,10) = 0; 
                        cha(zi,9) = cha(zi,9) + 1;
                        
                         % connecting to cell leader with correct weihgt
                        A(zi,zs) = 201;
                        A(zs, zi) = 201;
                        if(length(A) ~= n)
                         msg = 'Error occurred!!! Annoying 784';
                         error(msg)
                        end
                    else
                        
                        z3 = intersect(find(A(zi,:)==50), find(cha(:,1) == 2));
                        zp = length(z3);
                        if(zp > 5)
                            cha(zs,1) = 4;
                            cha(zs,10) = 0; 
                            cha(zi,9) = cha(zi,9) + 1;
                            % connect to other cell leader
                            A(zi,zs) = 201;
                            A(zs, zi) = 201;
                            
                             % Dividing cell
                            z4 = randperm(zp, 3);
                            z4 = z3(z4);
                            A(zi,z4) = 0;
                            A(z4, zi) = 0;
                            A(z4,zs) = 200;
                            A(zs, z4) = 50;
                            
                            z5 = setdiff(z3, z4);
                            A(z4, z5 ) = 0;
                            A(z5, z4 ) = 0;
                            A(z4, z4) = 50;
                            A(z5, z5) = 50;
                            
                        if(length(A) ~= n)
                         msg = 'Error occurred!!! Annoying 830';
                         error(msg)
                        end
                          
                            
                        else
                            cha(zs, 10) = cha(zs, 10) + rand/5;  
                        end
                        
                    end
                    
                else
                cha(zs,10) =  cha(zs, 10) -rand/5;
                cha(zi,9) = cha(zi,9) -1;
                
                        
                end
            else
             cha(zs, 10) = cha(zs, 10) + cha(zs,1)*(u/15)*(-1)^(cha(zs, 10) < u);  
                if(length(A) ~= n)
                msg = 'Error occurred!!! Annoying 353';
                    error(msg)
                end
            end
            %%%
            %%%
%%%%
%%%%%
%%%
            %L meets FS
         elseif( cha(zs,1) == 2)
             
             % Starting event
             if( A(zi,zs) == 200 && cha(zi,12) > 0 && cha(zs,13) < 1 )
                 even = even + 1;
                 z = [zi, find( A(zi,:)==200)];
                 
                % Commiting a crime is increased 
                cha(z, 5) = cha(z, 5) + 1;
                
                
                 
                 % This is finding military people connected to terrorist
                 [~, zm] = find( A(z,:));
                 zm = unique(zm);
                 zm = intersect(cha(:,1) ==1 , zm);
                 
                 y2 = [wT(1), wT(2).*ones(1, length(z)), wT(3).*ones(1, length(zm))]*[ cha(zi,12); cha(z, 9); -cha(zm,9)];
                 
                 %this is a probability distribution,which decides whether
                 %it is sucessful of not
                 y3 = (1/(1+exp((-y2 + muT)/sT)));
                 if(length(A) ~= n)
                  msg = 'Error occurred!!! Annoying 886';
                  error(msg)
                 end
                
  %%%%%%%%%%
  %%%%%%%%%%%
                 % This is the event is successful
                 if(y3>0.5)
                  suc = suc + 1;   
                 
                   
               % number of kills
                y = poissrnd(1) + 1;
                z3 = find(cha(:,1) == 0 | cha(:,1)==1);
                 z2 = randi( length(z3), 3,1 );
                 
                 z2 = z3(z2);
                nkill  = nkill +y;
                if(length(z2) >= y)
                    
                % increase exposure of violence
                 [ ~, b] = find(A(z2(1:y),:) > 0);
                    b = unique(b);
                    cha(b,6)  = cha(b,6) + 1; 
                    b = intersect(b, find(cha(:,1) ==0));
                           
                    % changing bias increase
                    cha(b,10 ) = cha(b,10 ) + .1 ;
                    
                    % Removing dead and repopulating resetpop(charact, x, F, edg, pop)
                        [ cha, A] = resetpop(cha, z2(1:y), A, n3, n);
                        z2 = z2(1:y);
                        
                else
                    
                 zl2 = length(z2);   
                y= y-zl2;
                for j2 = 1:zl2
                if(y>0)
                    z3 = setdiff(find(A(z2(j2),:)> 0),[ z'; z2]);
                  
                     if(~isempty(z3))
                         if(length(z3) >= y)
                                  
                % increase exposure of violence
                 [ ~, b] = find(A(z3(1:y),:) > 0);
                    b = unique(b);
                    cha(b,6)  = cha(b,6) + 1; 
                    b = intersect(b, find(cha(:,1) ==0));
                           if(length(A) ~= n)
                           msg = 'Error occurred!!! Annoying 945';
                             error(msg)
                            end
                    % changing bias increase
                    cha(b,10 ) = cha(b,10 ) + .1 ;
                        [ cha, A] = resetpop(cha, z3(1:y), A, n3, n);
                       
                         else
                            
                                  
                % increase exposure of violence
                 [ ~, b] = find(A(z3,:) > 0);
                    b = unique(b);
                    cha(b,6)  = cha(b,6) + 1; 
                    b = intersect(b, find(cha(:,1) ==0));
                    
                    % changing bias increase
                    cha(b,10 ) = cha(b,10 ) + .1 ;
                        [ cha, A] = resetpop(cha, z3, A, n3, n);
                       
                         end
                         y= y-length(z3);
                        
                     end
                     
                 
                  
                end
                  
                end 
                end
   
   
              % Terrorist reset
              cha(z, [9, 12]) = 0;
               cha(z,13) =  1; 
%%%%%%%
%%%%%%%%%%%   unsucesssful
                 else
                      % Lose all power if failed attempt
                    cha(z,[9,12]) = 0;
                    % wait time
                    cha(z,13) = 1;
                    
                     
                    % from UMD data 15.4 % of all unsucessful events ended with
                    % terrorist dieing
                    
                    if( rand < 0.155)
                        
                        
                    % increase exposure of violence
                 [ ~, b] = find( A(z,:) > 0 );
                    b = unique(b);
                    cha(b,6)  = cha(b,6) + 0.5; 
                        
                        % max number of terrorist that die
                       y = round(gamrnd(1,1)) + 1; 
                       
                       
                       zl = intersect(z,find(cha(:,1) == 4));
                       z = setdiff(z, zl);
                      if( y < length(z) )
                       tkill = tkill + y; 
                       
                       % Removing terrorist and repopulating resetpop(charact, x, F, edg, pop)
                        [ cha, A] = resetpop(cha, z(1:y), A, n3, n);
                            
                      else
                         
                          [ cha, A] = resetpop(cha, z, A, n3, n);
                          tkill = tkill + length(z) ;
                          if(length(A) ~= n)
                         msg = 'Error occurred!!! Annoying 1005';
                         error(msg)
                          end
                      end 
                    
                      
                    end
     
                 end
  %%%%
  %%%%%         end of event
  %%%%%
                
             else
          cha(zi,9) = cha(zi,9) + 1;   
          cha(zs,9) = cha(zs,9) + 1;  
          
             end
%%%
%%% Leader meets Military
         elseif( cha(zs,1) == 1 && (cha(zs,9)  > cha(zi,9) + 2))
             
            cha(zi, 13) =  1;
%%%
%%%
%%%         Leader meets Finacier   
         else
          cha(zi,9) = cha(zi,9) + 1;   
          cha(zs,9) = cha(zs,9) + 1;
         end
         
        end
 
 %%%%
 %%%%%
 %%%%
 %%%%%%%%%%%%%%%%%%%%%%%%
 %%%%
 %%%%
 %%%%%
 %%%%
 %%%%%%%%%%%%%%%%%%%%%%%%
 %%%%       
                  if(length(A) ~= n)
                    msg = 'Error occurred!!! Annoying 1017';
                   error(msg)
                  end
                  
         % Military 
    elseif(cha(zi,1) == 1)
        
        %meets Civilian
     if(cha(zs,1) == 0)
            % Increase or decrease bias towards military
        cha( zi , 11) = cha(zi, 11) + (u/5)*(-1)^(cha(zi, 11) < u);
        
        % Meets terrorist
     elseif(cha(zs,1) ~= 1)
         if(cha(zs,12) < arrest )
          cha(zs,9) = cha(zs,9) - 1;  
          cha(zi,9) = cha(zi,9) + 1;
          cha(zs, 13) = 1;
                       
          
         else
             % Total number captured
            cap = cap +1;
            
          cha(zi,9) = cha(zi,9) + 1;
             if( cha(zs,1) == 4)
               if( cha(zs,9) < cha(zi,9) +2)
                   % switch leader
                  a =  find(A(zs,:) == 200 );
                  for i2 = 1:length(a)
                a2(i2) =   1./(1+exp(-(wL*[cha(a(i2),2:7)'; (cha(a(i2),8) ==1)] - muL)/sL));
                  end 
                [ ~ ,b] = max(a2); 
                  cha(a(b),1) = 4;
                  [ cha, A] = resetpop(cha, zs, A, n3, n);
               else
              cha(zs,[9, 12]) = 0;
          cha(zs, 13) = 1;
               end
             elseif( cha(zs,9) < cha(zi,9)  )
                 
                  [ cha, A] = resetpop(cha, zs, A, n3, n);
             else 
                 
              cha(zs,[9, 12]) = 0;
              cha(zs, 13) = 1;
             end
                           
         end
         
         %M meets M
     else
         cha(zs,9) = cha(zs,9) + 1;
         cha(zi,9) = cha(zi,9) + 1;
                         
     end 
        
    end
                  
    end
end


% % % % Recruiting new military
% % if(mod(jj, 1000)==0)
% %    y = find(cha(:, 11) < 0 & cha(:, 10) > .9, 1);
% %    
% %    
% %     cha(y,1) = 1;
% %     cha(y,11) = 0;
% %     
% % end


%eachE(jj) = nkill;
end


ttday(k ,:) = [sum(cha(:,1)==2 ), sum(cha(:,1)==3),  sum(cha(:,1)==4 ),even, nkill, sum(cha(:,1)==1 ), tkill];

%tdayweight(:,:,10*(jjj-1)+k) = ttday;

end
toc
big(:,:, jjj) = ttday;
end
[15, 31, 35, ]
%%

% If you would like to find the averages for each thing shown
avB = zeros(40,6);
for i = 1:40
    
   avB(i,:) = sum(big(:,:,i))./50; 
    
    
end


