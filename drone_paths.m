clear all;
close all;
% starting from one point and scanning the area with minimum rotation
% fitness1: the most area to visit / the number of cells it can visit
% We want the value of the f1 function to be high

% fitness3: amount of change of direction
% We do not want the drones to make sharp turns.
% We want the value of the f3 function to be low 

% fitness2 : its proximity to the starting point
% Each drone needs to return to its starting point or as close as possible.


% Starting with P random solutions
% G times:
% generate P new solutions by crossover and mutation from solutions
% Crossover is made from 2 or 1 points according to the entered crossover value.
% mutation % mu
% Calculate the goodness values of the solutions produced
% Add the angles of rotation, add the number of different points traveled

% Calculate the distances of the generated solutions to the starting point and increase the probability of selecting values that approach the starting point.
% New generation: some of the old generation + some of the new generation


% possible moves [0 8] 
% 0 stand still
% 1-8 directions
hareket=[0 0; 0 1; -1 1; -1 0; -1 -1; 0 -1; 1 -1 ; 1 0; 1 1];

% cost between consecutive moves:
ac_t=[2:9];
ac=(ac_t-2)*45; ac=[0 ac]; % 0 0 45 90 135 180 225 270 315

G=200; % number of gen
P=100; % ppopulation size
cross=2; % crossover 1: single point crossover 2: from two point
bas=[1 1]; % starting point
drone_sayisi = 2; % number of drones

%  User will enter population size, gen value, crossover type, starting point and number of drones
prompt1 = 'Enter the population value';
prompt2 = 'Enter the generation value';
prompt3 = 'Enter the crossover type : 1 or 2';
prompt4 = 'Starting point : ';
prompt5 = 'Enter the number of total drone(s)';
P = input(prompt1);
G = input(prompt2);
cross = input(prompt3);
bas= input(prompt4);
drone_sayisi= input(prompt5);

Mx=9;My=Mx; % size of matrix 
M=zeros(Mx,My);% matrix
u=Mx*Mx-1; % len of each person 
mu=0.01; % mutation rate 
toplam_hamle = int32(P)/int32(drone_sayisi); % dividing the population for the number of drones to be used
toplam_hamle = toplam_hamle-1;
BK= toplam_hamle - toplam_hamle/2; % the best number of individuals directly copied to the next generation
B=zeros(P,u); % bireyler

son=[9 9]; % bitis noktasi
max_f1=u+1; % The maximum value of the fitness1 function determined to be used in normalization
max_f2=sum(abs(bas-son)); % The maximum value of the fitness2 function determined to be used in normalization
max_f3=180*(u-1);% The maximum value of the fitness3 function determined to be used in normalization

for z=1:drone_sayisi
    bb(z).best=zeros(1,G); % holds the degree of well-being of the best individual of each generation
    bb_yol(z).path=zeros(G,u); % keeps the path of the best individuals of every generation
end
ob=zeros(1,G); % keeps the average value of each generation
j=1;

artis = 0;
tic
% create first persons 
B=round(8*rand(P,u)); % fill with numbers between 0 and 8 

for i=1:G  % Values are used by generating new populations as much as the number of generations.
        % calculate the fitness values of individuals
    f1=zeros(1,toplam_hamle);
    for z=1:drone_sayisi   
        f3(z).value=f1(1:toplam_hamle); 
        f2(z).value=f1(1:toplam_hamle);
    end
    
    for j=1:toplam_hamle % bringing movements from the population
        M=zeros(Mx,My); % reset environment
        for z=1:drone_sayisi
            birey(z).array=B(j+artis,:);
            artis = artis + toplam_hamle;
        end
        artis = 0;
        % create the movement of the individual
        for a=1:drone_sayisi
            kxy(a).indis=bas;
        end
        M(kxy(1).indis,kxy(1).indis)=1;  % at start 
        
        for k=1:u % calculate movements of person 
            for z=1:drone_sayisi
                p_kxy(z).indis=kxy(z).indis+hareket(birey(z).array(k)+1,:); % make move 
            if p_kxy(z).indis(1)<1 || p_kxy(z).indis(2)<1 || p_kxy(z).indis(1)>Mx  || p_kxy(z).indis(2)>Mx % stand same if it goes out 
            else % didn't go out, move
                kxy(z)=p_kxy(z);
            end
            M(kxy(z).indis(1),kxy(z).indis(2))=k+1; % mark the point in the map 
            end
        end
        f1(j)=u-size(find(M),1)+1; % max number of cells visited
        for z=1:drone_sayisi
            f2(z).value(j)=sqrt(sum(abs(bas - kxy(z).indis) .^ 2)); % The proximity values to the starting point are calculated for each drone.
        %Calculation of deflection angle for each drone
            birey_arti=(birey(z).array)+1;
            birey_2=birey_arti(2:end);
            birey_1=birey_arti(1:end-1);
            a1=ac(birey_1); 
            a2=ac(birey_2);
            fark=abs(a1-a2);
            fark(fark>180)=360-fark(fark>180);
            f3(z).value(j)=sum(fark);
        end
    end
    % normalization for fitness functions 
        n_f1=f1/max_f1;
        for z=1:drone_sayisi
            n_f3(z).value=f3(z).value/max_f3;
            n_f2(z).value=f2(z).value/max_f2;
        end
        
        for z=1:drone_sayisi
            w(z).value =n_f1 + n_f2(z).value + n_f3(z).value; % sum of fitness values
            n_w(z).value=w(z).value/sum(w(z).value); % Mins must have a high chance of being selected in n_w
            yollar = B(:,:);
            n_w(z).value =1-n_w(z).value; % highlighting the most wanted cases that we want to be selected from the functions by inverting
            n_w(z).value=n_w(z).value/sum(n_w(z).value);
            % roulette wheel
			
            [sorted,inds]=sort(n_w(z).value);
            rn_w(z).value(inds)=1:toplam_hamle;
            rn_w(z).value=double(rn_w(z).value)/double(sum(rn_w(z).value));
            [val best_ind]=max(rn_w(z).value);
            %best_ind
            bb(z).best(i)=w(z).value(best_ind); % the best value is chosen
            bb_yol(z).path(i,:)=B(best_ind,:); % the path of the one with the best value is chosen
            ob(i)=mean(w(z).value);  
            secilenler(z).value = randsample(toplam_hamle,toplam_hamle,true,rn_w(z).value);% According to the values of the indices, indices with a high fitness function result are selected so that the chance of selection is higher.
            index_ayarla = (z-1) * toplam_hamle;
            secilenler(z).value = (secilenler(z).value) + double(index_ayarla);
        end
        % single and two point crossover according to the value entered to breed new individuals
        YB=zeros(P,u); % new individuals are created
        for p=1:drone_sayisi
            for j=1:toplam_hamle/2
                x = ((p-1) * toplam_hamle);
                b1=B(secilenler(p).value(j),:);
                b2=B(secilenler(p).value(j+(toplam_hamle/2)-1),:);
                if cross==1 % single point crossover
                    kesme=round((u-3)*rand(1,1))+2; % 2 - (u-1) between
                    YB(j+x,:)      =[b1(1:kesme) b2(kesme+1:end)];
                    YB(j+x+(toplam_hamle/2),:)=[b2(1:kesme) b1(kesme+1:end)];
                else
                     % two points crossover
                    kesme=round((u-3)*rand(1,2))+2; % 2 points between 2 - (u-1) 
                    kesme=sort(kesme); % sort in ascending order 
                    YB(j+x,:)=[b1(1:kesme(1)) b2(kesme(1)+1:kesme(2)) b1(kesme(2)+1:end) ];
                    YB(j+x+(toplam_hamle/2),:)=[b2(1:kesme(1)) b1(kesme(1)+1:kesme(2)) b2(kesme(2)+1:end) ];
                end
            end
        end
        if BK>0 % Copy the best BK value in B to YB
            YB(inds(BK+1:end),:)=B(inds(BK+1:end),:);
        end
        % apply mutation
        d_ind=rand(P,u)<mu; % cells to change
        yy=round(8*rand(P,u)); % what will they change
        YB(d_ind)=yy(d_ind);
        if BK>0 % Copy the best BK value in B to YB
            YB(inds(BK+1:end),:)=B(inds(BK+1:end),:);
        end
        
        B=YB; % new generation ready
        
    end

for a=1:drone_sayisi
    plot(bb(a).best); % draw the best
    hold on;
end
xlabel('Jenerasyon Sayisi');
ylabel('Fitness Degeri');
hold on 
plot(ob,'-r'); % draw the mean

for a=1:drone_sayisi % create paths for drawing the best selected drone movements
    [v bb_best(a).max] = max(bb(a).best);
    birey(a).yol = bb_yol(a).path(bb_best(a).max,:);
end
  M=zeros(Mx,My);
  for z=1:drone_sayisi
    kxy(z).indis =bas;
    rota(z).indis =zeros(2,u);
  end
  
  M(kxy(1).indis,kxy(1).indis)=1; % at start 

  for k=1:u  % Paths are created for the best selected drones
      for z=1:drone_sayisi
        p_kxy(z).indis=kxy(z).indis+hareket(birey(z).yol(k)+1,:);% create movement
      if  p_kxy(z).indis(1)<1 || p_kxy(z).indis(2)<1 || p_kxy(z).indis(1)>Mx  || p_kxy(z).indis(2)>Mx 
      else
        kxy(z)=p_kxy(z);
      end
        M(kxy(z).indis(1),kxy(z).indis(2))=k+1;
        rota(z).indis(:,k)=kxy(z).indis;
      end
end

%print the best individual
bb(end) % son en bireyin degeri
%M(son(1),son(2))=64;
% path olustu
figure;image((81/u)*M);
figure; plot(rota(1).indis(1,:),rota(1).indis(2,:));
hold on
for a=2:drone_sayisi
    plot(rota(a).indis(1,:),rota(a).indis(2,:));
    hold on
end


axis([0 Mx+1 0 Mx+1]);
toc      
