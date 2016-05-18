
clear all;
clc;
diary('mydata');
global N node Noise_Floor pktsTransmitted pktsLost;
N = 10;
%energy parameters
% Supply Voltage (V)
V=3;
% Listening Mode Current (Il)
Il=18.8*0.001;
% Transmit Mode Current (It)
It=17.4*0.001;
% Receive Mode Current (Ir)
Ir=19.7*0.001;
% Clear Channel Assessment (TCCA)
TCCA = 128*0.000001;
TCC = 128*0.00000001;
% Active Radio Current Iactive (Idle Listeing)
Ia=19.7*0.001;
% Byte Transmission Time (TB)
TB=32*0.000001;
% Sleep Current Is
%Switching Energy (Eidleswitch)
SE=827*0.000000001;
% Listen Mode Current (Ir)
Is= 0.2* 0.001;
%**************************************************************************
%Data Packet Length in bytes
P = 12;
% ACK Packet Length in Bytes
Ack = 11;
%number of packts
m = 10;
%**************************************************************************
%Time for one symbol in micro-seconds
ST = 0.000016;
%Byte Transmission Time in microseconds
TB = 32*0.000001;
%Backoff period in microseconds
BP = 0.00032;
% Clear Channel Assessment Time(TCCA) in microsec
TCCA = 128*0.000001;
%The minimum value of the backoff exponent (BE) in the CSMA-CA algorithm
macMinBE=3;
%The maximum value of the backoff exponent (BE) in the CSMA-CA algorithm
macMaxBE=5;
%The maximum number of backoffs the CSMA-CA algorithm will attempt before
%declaring a channel access failure
macMaxCSMABackoffs=4;
% RX-to-TX or TX-to-RX maximum turnaround time (in microseconds0
aTurnaroudTime = 0.000192;
%CAP period start time
CAPStartTime = 0;
%CAP period end time
CAPEndTime = 15360; %in microseconds
%SIFS
SIFS = 0.000192;
% The maximum time to wait for an acknowledgment frame to
% arrive following a transmitted data frame in microseconds
macACKWaitDuration = 0.00192;
MIN = 2;%default in zigbee 2
Noise_Floor = -100;
%**************************************************************************
CW = [2 2 2 2 2 2 2 2 2 2];
NB = zeros(1,N);
backoff = zeros(1,N);
    Isources = [];
pktsTransmitted = zeros(1,N);
pktsLost = zeros(1,N);
E1 = zeros(1,N);
pktsTransmitted(1,1) = 1;
N_point = 10;
D = zeros(1,N);
node = [0.5 1.2; rand 2*rand ;rand 2*rand ; rand 2*rand ;rand 2*rand ;rand 2*rand ;rand 2*rand ;rand 2*rand ;rand 2*rand ;rand 2*rand];
for j=1:N
    D(1,j) =  topo_dist(1, j);
end
counter = 1;
while counter < 10
    NB = zeros(1,N);
    csmaTime = zeros(1,N);
    pktStartTXTime = zeros(1,N);
    pktEndTXTime = zeros(1,N);
    baseChannel = zeros(1,N);
    
    sourceState = floor(2*rand(1,N));
    sourceState(1) = 20;
    
    sources = find(sourceState>=1 & sourceState < 20);
    S = length(sources); %sources size
    
    backoff = floor(2^MIN*rand(1,S));
    
    for s=1:S
        Isources(s).srcID = sources(s);
        i = 1;
        for t =1:S
            if backoff(s) == backoff(t)
                Isources(s).b(i) = sources(t);
                i = i + 1;
            end
        end
    end
    %**************************************************************************
    %create an array structure of each source and its corresponding
    %interferers
    tr = zeros(1,S);
    for s=1:S
        if length(Isources(s).b) == 1
            tr(s) = sources(s);
        end
    end
    %scheduling the event by ascending sorting the backoff array
    sortMatrix = zeros(S,3);
    for i = 1:S
        sortMatrix(i,1) = backoff(i);
        sortMatrix(i,2) = sources(i);
        sortMatrix(i,3) = D(sources(i))*100;
    end
    sortMatrix = sortrows(sortMatrix,1);
    for i = 1:S
        backoff(i)= sortMatrix(i,1);
        sources(i)= sortMatrix(i,2);
    end
    fprintf('\nSources = ');
    disp(sources);
    fprintf('\nBackoff = ');
    disp(backoff);
    src2Cdist = sortMatrix(:,3)';
    %**************************************************************************
    %number of transmitters on basechannel;
    trans = find(tr>0);
    sizeT = length(trans);
    if sizeT == 0
        continue;
    end
    transmitters = zeros(1,sizeT);
    sortedTransmitters = zeros(1,sizeT);
    sortedBackoff = zeros(1,sizeT);
    p = 1;
    for s=1:S
        if tr(s) > 0
            transmitters(p) = tr(s);
            p = p + 1;
        end
    end
    [c,ia]=intersect(sources,transmitters);
    trans = sort(ia);
    sortedTransmitters = sources(trans);
    fprintf('Sorted transmitters are: ');
    disp(sortedTransmitters);
    fprintf('Sorted backoff: ');
    sortedBackoff = backoff(trans);
    disp(sortedBackoff);
    sizeSB = length(sortedBackoff);
    maxBO = sortedBackoff(sizeSB);
    %for all transmitters set the channel as clear
    for  f = 1:sizeT
        baseChannel(sortedTransmitters(f)) = 1;
    end
    fprintf('\nBase channel : ');
    disp( baseChannel);

    %initialize the sources backoff array
    for g =1:S
        csmaTime(g) = backoff(g);
    end
    %*********************************************************************

    %Initialize sensors energy
    for i = 1:N
        if i == 1 % coordinator
            E1(i)= 0;
        else
            E1(i)=E1(i)+(2*SE);
        end
    end
    for i = 1:S
        src = sources(i);
        while NB(src) <= macMaxCSMABackoffs
            csmaTime(src) = csmaTime(src) + 1;
            if baseChannel(src)
                CW(src) = CW(src) - 1;
                if CW(i) == 0
                    fprintf('Source %d sends data packet \n',src);
                else
                    csmaTime(src) = csmaTime(src) + 1;
                end
            else
                CW(src) = 2;
                NB(src) = NB(src) + 1;
                backoff(src) = min(backoff(src) + 1, macMaxBE);
                if NB(src) <= macMaxCSMABackoffs
                    backoff(src) = floor(2^backoff(i)*rand());
                    csmaTime(src) = backoff(i) + 1;
                else
                    % Algorithm terminates with failure
                end
            end

        end
        counter = counter + 1;
        if counter == 10
            break;
        end
    end
end


