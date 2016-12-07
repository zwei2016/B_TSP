% classic TSP
% short range aircraft + european cities 
% Let's imagie one senario : Airbus sends its Beluga aircraft (A300-600ST)
% to carry the oversized cargo from one european airport to Toulouse. 
% As a heavy transporter, the range of Beluga is only 2779 km.
% the trips for the 40 airports is below 9000 kw in general 


clear all;

% 40 airports selected from the list of the busiest airports in Europe. 
airport_name = {
    'Toulouse';    %1
    'Paris';
    'Lyon';        %3
    'Brussels';
    'Prague';      %5
    'Copenhagen';
    'Helsinki';
    'Frankfurt';
    'Dusseldorf';
    'Berlin';      %10
    'Thessaloniki';
    'Larnaca';
    'Budapest';
    'Reykjavik';
    'Dublin';      %15
    'Naples';
    'Milan';
    'Oslo';
    'Warsow';
    'Lisbon';      %20
    'Moscow';
    'Belgrade';
    'Madrid';      %23
    'Barcelona';   %24
    'Seville';     %25
    'Geneva';
    'Zurich';
    'Istanbul';
    'Izmir';
    'London';      %30
    'Amsterdam';
    'Ankara';
    'Pisa';
    'Belfast';     %34
    'Vienna';
    'Nantes';      %36
    'Munich';
    'Hamburg';     %38
    'Marseille';
    'Bologna';     %40
    };

N = size(airport_name);

% check if the airport position file exits
% if OK input the info from that file
% if no OK just creat the file with python 
fid = fopen('position_info.txt','r');
if fid <0
    % output the name list in airport_list.txt
    file_airport_list = fopen('airport_list.txt','w');
    for i=1:N(1)
        fprintf(file_airport_list,'%s \n', airport_name{i});
    end
    fclose(file_airport_list);
    % call the python file: geolocation.py, using Python Geocoding Toolbox geopy 1.11.0
    %[status, cmdout] = system('python geolocation.py airport_list.txt position_info.txt')
    [status, cmdout] = system('python geolocation.py airport_list.txt position_info.txt');
    pause(15);
end

formatSpec = '%s %f %f';
% input the airport position from position_info.txt
airport = struct ('city', [], 'latitude',[], 'longitude', []);
[c,lat,lon] = textread('position_info.txt',formatSpec);
fclose('all');
    
for i=1:N(1)
    airport(i).city = c(i);
    airport(i).latitude = lat(i);
    airport(i).longitude = lon(i);
end

% after the preparation, we can calcalat the distance between every two
% cities (from city i to city j), in order to the next optimization problem

Dis_cities = zeros(N(1));
for i=1:N(1)
    for j=1:N(1)
        if i==j
            Dis_cities(i,j) = 0;
        else
            Dis_cities(i,j) = getdistance(airport(i).latitude,airport(i).longitude,airport(j).latitude,airport(j).longitude);
        end
    end
end

figure(1)
% draw the european map with the cities 
h = worldmap('Europe');
%getm(h,'MapProjection');
getm(h,'MapParallels');
geoshow('landareas.shp', 'FaceColor', [0.15 0.5 0.15])
geoshow('worldlakes.shp', 'FaceColor', 'cyan')
%geoshow('worldrivers.shp', 'Color', 'blue')
%geoshow('worldcities.shp', 'Marker', '.', 'Color', 'red')
geoshow(lat,lon,'DisplayType','point')

for i=1:N(1)
   textm(airport(i).latitude, airport(i).longitude, airport(i).city) 
end

% Map is OK, now is the optimization problem
nStops = N(1);
% Generate all the trips, meaning all pairs of stops.
idxs = nchoosek(1:nStops,2);

%dist = getdistance (lat(idxs(:,1)),lon(idxs(:,1)),lat(idxs(:,2)),lon(idxs(:,2)));
%dist = hypot(lat(idxs(:,1)) - lat(idxs(:,2)), lon(idxs(:,1)) - lon(idxs(:,2)));
account = 1;
dist = [];
for i = 1:nStops 
    for j = i+1: nStops
        % change the matrix to list
        dist(account) = Dis_cities(i,j);
        account = account + 1;
    end
end

lendist = length(dist);

Aeq = spones(1:length(idxs));
beq = nStops;

Aeq = [Aeq;spalloc(nStops,length(idxs),nStops*(nStops-1))]; % allocate a sparse matrix

for ii = 1:nStops
    whichIdxs = (idxs == ii); % find the trips that include stop ii
    whichIdxs = sparse(sum(whichIdxs,2)); % include trips where ii is at either end
    Aeq(ii+1,:) = whichIdxs'; % include in the constraint matrix
end
beq = [beq; 2*ones(nStops,1)];

intcon = 1:lendist;
lb = zeros(lendist,1);
ub = ones(lendist,1);

opts = optimoptions('intlinprog','Display','off');
[x_tsp,costopt,exitflag,output] = intlinprog(dist,intcon,[],[],Aeq,beq,lb,ub,opts);

tours = detectSubtours(x_tsp,idxs);
numtours = length(tours); % number of subtours
fprintf('# of subtours: %d\n',numtours);

%%%%
A = spalloc(0,lendist,0); % Allocate a sparse linear inequality constraint matrix
b = [];
while numtours > 1 % repeat until there is just one subtour
    % Add the subtour constraints
    b = [b;zeros(numtours,1)]; % allocate b
    A = [A;spalloc(numtours,lendist,nStops)]; % a guess at how many nonzeros to allocate
    for ii = 1:numtours
        rowIdx = size(A,1)+1; % Counter for indexing
        subTourIdx = tours{ii}; % Extract the current subtour
%         The next lines find all of the variables associated with the
%         particular subtour, then add an inequality constraint to prohibit
%         that subtour and all subtours that use those stops.
        variations = nchoosek(1:length(subTourIdx),2);
        for jj = 1:length(variations)
            whichVar = (sum(idxs==subTourIdx(variations(jj,1)),2)) & ...
                       (sum(idxs==subTourIdx(variations(jj,2)),2));
            A(rowIdx,whichVar) = 1;
        end
        b(rowIdx) = length(subTourIdx)-1; % One less trip than subtour stops
    end

    % Try to optimize again
    [x_tsp,costopt,exitflag,output] = intlinprog(dist,intcon,A,b,Aeq,beq,lb,ub,opts);
        % How many subtours this time?
    tours = detectSubtours(x_tsp,idxs);
    numtours = length(tours); % number of subtours
    fprintf('# of subtours: %d\n',numtours);
    %segments = find(x_tsp);
    %lh = zeros(nStops,1);
    %lh = updateSalesmanPlot(lh,x_tsp,idxs,lon,lat);
    %geoshow(lat,lon,x_tsp,'DisplayType','point')
    
end

%%%%%%%%%%%%%
% Update the map with the Beluga iteration 
%%%%%%%%%%%%%%
Res = zeros(nStops);
account = 0;
for i = 1: nStops 
   for j = i+1: nStops
     account = account + 1;
     if x_tsp(account)==1;
        Res(i,j) = 1;

        [pointlat,pointlon] = gcwaypts(lat(i),lon(i),lat(j),lon(j));
        geoshow(pointlat,pointlon,'DisplayType','line')
        pointTrack = track2('gc',lat(i),lon(i),lat(j),lon(j));
        plotm(pointTrack);
     end 
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Biograph of the iteration 
%%%%%%%%%%%%%%%%
ListA = [];
ListB = [];
ListW = [];
Edge = struct ('cityA', [], 'cityB',[]);
flag = 0;
for i = 1: nStops
  for j = 1: nStops
    if Res(i,j)==1
       flag=flag+1;
       ListA(flag)=i;
       ListB(flag)=j;
       ListW(flag)= Dis_cities(i,j);
       Edge(flag).cityA = c(i) ;               % edge connecting two cities
       Edge(flag).cityB = c(j);
    end
  end
end

% several graph to illustrate the iteration 
DG = sparse(ListA,ListB,ListW,nStops,nStops);
UG = tril(DG + DG');
Graph = biograph(UG,airport_name,'ShowArrows','off','ShowWeights','on')
view(Graph);
%isdag(Graph);

% using DFS to find the circle 
figure(2)
G = graph(ListA,ListB);
plot(G);
title('the iteration of aircraft');
% get the list of iternation of aircraft
Iter = dfsearch(G,1);
% Generate the Markov Chain Transition Matrix from Iter
% A node of G is a state of MDP
% An arc of G is an action in MDP: deterministic arc
TransMDP = zeros(nStops)
for i = 1:nStops-1
    TransMDP(Iter(i),Iter(i+1))=1;
end
TransMDP(Iter(nStops),1)=1;
       
% using the hamiltonian to get the hamilton path and the visulized path
% basing on the MDP transition matrix
% idea is from Hamiltonian cycle problem (HCP)
figure(3)
hamPath = hamiltonian(TransMDP,1,1);
ListS = hamPath(1:length(hamPath)-1);
ListT = hamPath(2:length(hamPath));
EdgeTable = table([ListS' ListT'],'VariableNames',{'EndNodes'});
% creat the EdgeTable from the hamPath
plot(graph(EdgeTable));
title('Hamiltonian path for routine mission');

%%%%%%%%%%%%%%%%%
%%factory sites trip: Nantes 36, Hamburg 38, Madrid 23, Seville 25 
%%% two potential sites: Lyon 3,  Barcelona 24
%%% Beluga will probably often visit these sites directly from Toulouse 1.
%%% Therefore we add consider a new matrix Mission with M(1,36),M(36,1)
%%% M(1,38), M(38,1), M(1,23), M(23,1), M(1,25), M(25,1), M(1,3), M(3,1)
%%% M(1,24), M(24,1) are equal to 1.
%%%%%%%%%%%%%
M = zeros(nStops);
Fac_sites = find(ismember(airport_name,{
    'Nantes','Hamburg','Madrid','Seville','Lyon','Barcelona'}))
M(1,Fac_sites)=1;
M(Fac_sites,1)=1;
% The routine mission and the special one can be done in one trip as well
Total_M1 = double(TransMDP | M);
% This is not a MDP transition matris, but it is a graph
figure(4)
hamPath = hamiltonian(Total_M1,1,1);
ListS = hamPath(1:length(hamPath)-1);
ListT = hamPath(2:length(hamPath));
EdgeTable = table([ListS' ListT'],'VariableNames',{'EndNodes'});
plot(graph(EdgeTable));
title('Hamiltonian path for two missions');
%%%%normally the Hamiltonian path stays the same as previous

%%%%%%%%%
%%suppose the special mission only has one third frequence than the routine
%%one: Pr(s)=0.25 and Pr(r)=0.75
Total_M2 = 0.25*M + 0.75*TransMDP;

for i=1:N(1)
    if sum(Total_M2(i,:))~=1
        Total_M2(i,:)= normalization(Total_M2(i,:));
    end
end

Num = 150;
row = size(Total_M2,1);
x_states = zeros(row,Num);
x = eye(1,row);
x_states(:,1) = x;

% after Num 100 missions completed in general 
for t = 1:Num %col
    x = double(x*Total_M2);
    
    x_states(:,t+1) = x; 
end
figure(5);
plot(1:N(1),x_states(:,end),'-b','linewidth',4);
hold on
stem(1:N(1),x_states(:,end));
grid
xlabel('40 airports in Europe');
ylabel('Probability Distribution after 150 missions');
title('Probability and Importance of Airport');
set(gca,'fontsize',18)
axis tight

%% steady state
[V, l] = eigs(Total_M2',1,'lm');
%the columns of V are the eigenvectors
steady_state_value = normalization(V);
figure(6)
plot(1:40,steady_state_value,'-r','linewidth',3);
grid
xlabel('40 airports in Europe');
ylabel('Probability Distribution');
title('Steady State Probability Distribution Values');
set(gca,'fontsize',18)
axis tight









