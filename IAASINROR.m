
clear all;
clc;
%diary('mydata');
global N node relayState Noise_Floor pktsTransmitted pktsLost;
N = 11;

B = 20 * 1000; %Bandwidth in Hertz
symbol = 6; %number of bits in one symbol
L = 2^symbol;
capacity = 2 * B * log(L)/log(2); % maximum bit rate
acheivableBitRate = floor(capacity/2);
sinrThr = 2^((log(2)*acheivableBitRate) / B) - 1;
disp(sinrThr);
disp(acheivableBitRate);


pktsTransmitted = zeros(1,N);
RpktsTransmitted =  zeros(1,N);
pktsLost = zeros(1,N);
E1 = zeros(1,N);
E3 = zeros(1,N);
pktsTransmitted(1,1) = 1;
relayEnergy = zeros(1,N);
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
macMinBE=4;
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
r = 90;
Noise_Floor = -100;
%**************************************************************************
N_point = 10;
D = zeros(1,N);
%node = [0.5 1.2; rand 2*rand; rand 2*rand; rand 2*rand; rand 2*rand; rand 2*rand; rand 2*rand; rand 2*rand; rand 2*rand; rand 2*rand; rand 2*rand ; rand 2*rand; rand 2*rand; rand 2*rand; rand 2*rand; rand 2*rand ; rand 2*rand; rand 2*rand; rand 2*rand; rand 2*rand];%20 nodes
%node = [0.5 1.2; rand 2*rand; rand 2*rand; rand 2*rand; rand 2*rand; rand 2*rand; rand 2*rand; rand 2*rand; rand 2*rand; rand 2*rand; rand 2*rand ; rand 2*rand; rand 2*rand; rand 2*rand; rand 2*rand];%15 nodes
node = [0.5 1.2; 0.4 1.9; 0.2 0.5; 0.4 0.6; 0.3 1.4; 0.2 1.9; 0.6 0.3; 0.7 0.8; 0.55 0.9; 0.5 1.6; 0.6 1.9]; %10 nodes
% plot(node(:,1),node(:,2), 'b*')
% hold
% grid on
for j=1:N
    D(1,j) =  topo_dist(1, j);
end
counter = 1;
rc1 = zeros(1,N);
numberMessagesTransmitted = 0;
SINRthr = 20;
srcPktsTXed = zeros(1,N);

while counter < 1000
    


    %fprintf('\n/*** -- NEW CAP PERIOD -- ***/\n');
    relayState = zeros(1,N);
    interCaseRelays = zeros(1,N);
    successfulRelays = zeros(1,N);
    reservedCaseRelays = zeros(1,N);
    interRelays = zeros(1,N);
    succRelays = zeros(1,N);
    resRelays = zeros(1,N);
    Isources = [];
    NB = zeros(1,N);
    rc1 = zeros(1,N);
    csmaTime = zeros(1,N);
    pktStartTXTime = zeros(1,N);
    pktEndTXTime = zeros(1,N);
    baseChannel = zeros(1,N);
    CW = [2 2 2 2 2 2 2 2 2 2];
    % sourceState = floor(3*rand(1,N));
    sourceState(1) = 20;
    sourceState(2) = 1;
    sourceState(3) = 1;
    sourceState(4) = 0;
    sourceState(5) = 0;%relay
    sourceState(6) = 1;%source
    sourceState(7) = 1;
    sourceState(8) = 1;
    sourceState(9) = 0;
    sourceState(10) = 0;
    sourceState(11) = 1;
    minBE = 4;
    sources = find(sourceState>=1 & sourceState < 20);
    idleRelays = find(sourceState == 0);
    sizeR = length(idleRelays);
    S = length(sources); %sources size
    backoff = floor(2^minBE*rand(1,S));
    %for each source, find its interfering sources and save in a structure
    %array
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
    sizeR = length(idleRelays); %all idle relays size and coordinator is excluded
    fprintf('\nBackoff = ');
    disp(backoff);
    fprintf('\nIdle relays = ');
    disp(idleRelays);
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
    D = D * 100;
    sizeSB = length(sortedBackoff);
    maxBO = sortedBackoff(sizeSB);
    %for all transmitters set the channel as clear
    for  f = 1:sizeT
        baseChannel(sortedTransmitters(f)) = 1;
    end
    fprintf('\nBase channel : ');
    disp( baseChannel);
    %**********************************************************************
    intSrcs = setxor(sources,transmitters); %complement of transmitters in sources
    sizeIS = length(intSrcs);
    fprintf('\nInterfering sources are: ');
    disp(intSrcs);
    src2relayD = zeros(1,sizeR);
    ISBackoff  = floor(2^minBE*rand(1,sizeIS));
    ISBackoff = ISBackoff + maxBO;
    [c,indCols] = sort(ISBackoff);
    sortedInterferingSources = intSrcs(indCols);
    fprintf('\nSorted Interfering Sources:');
    disp(sortedInterferingSources);
    fprintf('\nSorted ISBackof:');
    disp(ISBackoff(indCols));
    ISBaseChannel = zeros(1,sizeIS);

    %**********************************************************************
    SINRthreshold = 100;
    sizeIS = S - sizeT; % size of interfering sources
    src2rlyDist = [];
    %array holds the powers received from a given transmitter at all relays
    P_Received = zeros(sizeT,sizeR);

    %initialize the sources backoff array
    for g =1:S
        csmaTime(g) = backoff(g);
    end
    %**************************************************************************
    %find distance from each relay to coordinator
    for k=1:sizeR % size of all relays
        rly2coorDist(k) = topo_dist(1,idleRelays(k)); % a source to relays distance
    end
    %**************************************************************************
    relayInterferences = zeros(sizeT,sizeR);
    relaysSINR = zeros(sizeT,sizeR);
    for i = 1:S % for all sources
        src = sources(i);
        d = topo_dist(sources(i),1)*100;
        lastTransmitter = sortedTransmitters(sizeSB);
        lastBackoff = maxBO;
        if ismember(src,transmitters) % checks if the source has clear channel to transmit (i.e a transmitter)
            if d <= r %the source is within a distance d from coordinator
                NB(src) = 0;
                CW(src) = 2;
                csmaTime(src) = backoff(i)*BP + TCCA;
                E1(src)=E1(src)+ TCCA*Il*V; % Energy consumption for whole WBAN
                E3(src)=E3(src)+ TCCA*Il*V; %  Energy consumption for sources only
                fprintf('\nBase Channel Case: Node %d performs CCA', src);
                while NB(src) <= macMaxCSMABackoffs
                    if baseChannel(src) == 1
                        CW(src) = CW(src) - 1;
                        if CW(src) == 0
                            numberMessagesTransmitted = numberMessagesTransmitted + 1;
                            baseChannel(src) = 0;
                            CW(src) = 2;
                            NB(src) = 0;
                            %Source transmits m data packets
                            pktStartTXTime(src) = csmaTime(src) + ST + SIFS; % pkts start transmission time
                            pktEndTXTime(src) = pktStartTXTime(src) + m*P*TB; % pkts end transmission time
                            fprintf('\nBase Channel Case: Node %d transmits data packet', src);
                            E1(src)=E1(src) + (m*P*TB*It*V); %energy cost for transmiting m packets
                            %E1(src)=E1(src) + (macACKWaitDuration*Il*V); % %energy cost to wait (listening) an ack for the last packet transmitted
                            %E1(src)=E1(src) + (Ack*TB*Ir*V); % energy cost for receiving one ack
                            %E3(src)=E3(src) + (m*P*TB*It*V); %energy cost for transmiting m packets
                            %E3(src)=E3(src) + (macACKWaitDuration*Il*V); % %energy cost to wait (listening) an ack for the last packet transmitted
                            %E3(src)=E3(src) + (Ack*TB*Ir*V); % energy cost for receiving one ack
                            pktsTransmitted(src) = pktsTransmitted(src) + 1;
                            srcPktsTXed(src) = srcPktsTXed(src) + 1;
                            %compute powers received at each of all relays
                            for l = 1:sizeT
                                for k=1:sizeR% size of all relays
                                    s2rdist = topo_dist(transmitters(l),idleRelays(k)); % a source to relays distance    %%% add SINR instead
                                    src2rlyDist(k) = s2rdist;
                                    [pow_received, Path_Loss, snr]  = Implant_PL_Model(s2rdist,N_point,2);
                                    P_Received(l,k) = pow_received;
                                    fprintf('\n');
                                    disp( P_Received(l,k));
                                end
                            end
                            %Compute all relays SINRs
                            maxRelaySINR = zeros(1,sizeT);
                            for l = 1:sizeT
                                for k =1:sizeR
                                    %relayInterferences(l,k) =  sum(P_Received(:,k)) - P_Received(l,k);
                                    relayInterferences(l,k) = 0; % this is because sources do not produce interference on the relays ince they are at least 1
                                    relaysSINR(l,k) = P_Received(l,k) - Noise_Floor; %- relayInterferences(l,k)%
                                end
                            end
                            fprintf('\nBase Case; Relays SINR;');
                            disp(relaysSINR(1,:));
                            avg = mean(relaysSINR(l,:));
                            s = std(relaysSINR(l,:));
                            b = find(relaysSINR(l,:) >= sinrThr);
                            [val, ind]= max(relaysSINR(l,:));
                            idx = idleRelays(ind);
                            if length(b) == 1
                                if ~ismember(idx,successfulRelays) && relayState(idx)== 0
                                    pktStartTXTime(idx) = pktEndTXTime(src) + ST + aTurnaroudTime + SIFS ; % pkts start transmission time
                                    pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkts end transmission time
                                    fprintf('\nBase Channel Case: Relay %d of Source %d replies with an Ack',idx,src);
                                    baseChannel(src) = 0;
                                    CW(src) = 2;
                                    NB(src) = 0;
                                    successfulRelays(idx) = idx;
                                    relayState(idx) = 0;
                                    E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                    %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                    relayEnergy(idx) = E1(idx);
                                    pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                    RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                    break;
                                else
                                    fprintf('\nBase Channel Case: Source %d could not find any relays, "Packet lost"', src);
                                    baseChannel(src) = 0;
                                    CW(src) = 2;
                                    NB(src) = 0;
                                    pktsLost(src) = pktsLost(src) + 1;
                                    break;
                                end
                            else
                                if length(b) == 2
                                    [val, ind]= max(relaysSINR(l,:));
                                    idx = idleRelays(ind);
                                    if ~ismember(idx,successfulRelays) && relayState(idx)== 0
                                        pktStartTXTime(idx) = pktEndTXTime(src) + ST + aTurnaroudTime + SIFS; % pkts start transmission time
                                        pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkts end transmission time
                                        fprintf('\nBase Channel Case: Relay %d of Source %d replies with an Ack',idx,src);
                                        baseChannel(src) = 0;
                                        CW(src) = 2;
                                        NB(src) = 0;
                                        relayState(idx) = 0;
                                        successfulRelays(idx) = idx;
                                        E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                        %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                        relayEnergy(idx) = E1(idx);
                                        pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                        RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                        break;
                                    else
                                        relaysSINR(l,ind) = -500;
                                        [secondMax, secondIndex] = max(relaysSINR(l,:));
                                        idx = idleRelays(secondIndex);
                                        if ~ismember(idx,successfulRelays) && relayState(idx)== 0
                                            pktStartTXTime(idx) = pktEndTXTime(src) + ST + aTurnaroudTime + SIFS; % pkts start transmission time
                                            pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkts end transmission time
                                            fprintf('\nBase Channel Case: Relay %d of Source %d replies with an Ack',idx,src);
                                            baseChannel(src) = 0;
                                            CW(src) = 2;
                                            NB(src) = 0;
                                            relayState(idx) = 0;
                                            successfulRelays(idx) = idx;
                                            E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                            %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                            relayEnergy(idx) = E1(idx);
                                            pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                            RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                            break;
                                        else
                                            fprintf('\nBase Channel Case: Source %d could not find any relays, "Packet lost"', src);
                                            baseChannel(src) = 0;
                                            CW(src) = 2;
                                            NB(src) = 0;
                                            pktsLost(src) = pktsLost(src) + 1;
                                            break;
                                        end
                                    end
                                else
                                    if length(b) == 3
                                        [val, ind] = max(relaysSINR(l,:));
                                        idx = idleRelays(ind);
                                        if ~ismember(idx,successfulRelays) && relayState(idx)== 0
                                            pktStartTXTime(idx) = pktEndTXTime(src) + ST + aTurnaroudTime + SIFS; % pkts start transmission time
                                            pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkts end transmission time
                                            fprintf('\nBase Channel Case: Relay %d of Source %d replies with an Ack',idx,src);
                                            baseChannel(src) = 0;
                                            successfulRelays(idx) = idx;
                                            CW(src) = 2;
                                            NB(src) = 0;
                                            relayState(idx) = 0;
                                            successfulRelays(idx) = idx;
                                            %E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                            E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                            relayEnergy(idx) = E1(idx);
                                            pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                            RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                            break;
                                        else
                                            relaysSINR(l,ind) = -500;
                                            [secondMax, secondIndex] = max(relaysSINR(l,:));
                                            idx = idleRelays(secondIndex);
                                            if ~ismember(idx,successfulRelays) && relayState(idx)== 0
                                                pktStartTXTime(idx) = pktEndTXTime(src) + ST + aTurnaroudTime + SIFS; % pkts start transmission time
                                                pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkts end transmission time
                                                fprintf('\nBase Channel Case: Relay %d of Source %d replies with an Ack',idx,src);
                                                baseChannel(src) = 0;
                                                CW(src) = 2;
                                                NB(src) = 0;
                                                relayState(idx) = 0;
                                                successfulRelays(idx) = idx;
                                                E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                relayEnergy(idx) = E1(idx);
                                                pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                break;
                                            else
                                                relaysSINR(l,secondIndex) = -500;
                                                [thirdMax, thirdIndex] = max(relaysSINR(l,:));
                                                idx = idleRelays(thirdIndex);
                                                if ~ismember(idx,successfulRelays) && relayState(idx)== 0
                                                    pktStartTXTime(idx) = pktEndTXTime(src) + ST + aTurnaroudTime + SIFS; % pkts start transmission time
                                                    pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkts end transmission time
                                                    fprintf('\nBase Channel Case: Relay %d of Source %d replies with an Ack',idx,src);
                                                    baseChannel(src) = 0;
                                                    CW(src) = 2;
                                                    NB(src) = 0;
                                                    relayState(idx) = 0;
                                                    successfulRelays(idx) = idx;
                                                    E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                    %E1(idx)=E1(idx) + Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                    relayEnergy(idx) = E1(idx);
                                                    pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                    RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                    break;
                                                else
                                                    fprintf('\nBase Channel Case: Source %d could not find any relays, "Packet lost"', src);
                                                    baseChannel(src) = 0;
                                                    CW(src) = 2;
                                                    NB(src) = 0;
                                                    pktsLost(src) = pktsLost(src) + 1;
                                                    break;
                                                end
                                            end
                                        end
                                    else
                                        if length(b) == 4
                                            [val, ind] = max(relaysSINR(l,:));
                                            idx = idleRelays(ind);
                                            if ~ismember(idx,successfulRelays) && relayState(idx)== 0
                                                pktStartTXTime(idx) = pktEndTXTime(src) + ST  + aTurnaroudTime + SIFS; % pkt start transmission time
                                                pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                                fprintf('\nBase Channel Case: Relay %d of Source %d replies with an Ack',idx,src);
                                                baseChannel(src) = 0;
                                                CW(src) = 2;
                                                NB(src) = 0;
                                                relayState(idx) = 0;
                                                successfulRelays(idx) = idx;
                                                E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                %E1(idx)=E1(idx) + Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                relayEnergy(idx) = E1(idx);
                                                pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                break;
                                            else
                                                relaysSINR(l,ind) = -500;
                                                [secondMax, secondIndex] = max(relaysSINR(l,:));
                                                idx = idleRelays(secondIndex);
                                                if ~ismember(idx,successfulRelays) && relayState(idx)== 0
                                                    pktStartTXTime(idx) = pktEndTXTime(src) + ST  + aTurnaroudTime + SIFS; % pkt start transmission time
                                                    pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                                    fprintf('\nBase Channel Case: Relay %d of Source %d replies with an Ack',idx,src);
                                                    baseChannel(src) = 0;
                                                    CW(src) = 2;
                                                    NB(src) = 0;
                                                    relayState(idx) = 0;
                                                    successfulRelays(idx) = idx;
                                                    E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                    %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                    relayEnergy(idx) = E1(idx);
                                                    pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                    RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                    break;
                                                else
                                                    relaysSINR(l,secondIndex) = -500;
                                                    [thirdMax, thirdIndex] = max(relaysSINR(l,:));
                                                    idx = idleRelays(thirdIndex);
                                                    if ~ismember(idx,successfulRelays) && relayState(idx)== 0
                                                        pktStartTXTime(idx) = pktEndTXTime(src) + ST  + aTurnaroudTime + SIFS; % pkt start transmission time
                                                        pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                                        fprintf('\nBase Channel Case: Relay %d of Source %d replies with an Ack',idx,src);
                                                        baseChannel(src) = 0;
                                                        CW(src) = 2;
                                                        NB(src) = 0;
                                                        relayState(idx) = 0;
                                                        successfulRelays(idx) = idx;
                                                        E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                        %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                        relayEnergy(idx) = E1(idx);
                                                        pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                        RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                        break;
                                                    else
                                                        relaysSINR(l,thirdIndex) = -500;
                                                        [fourthMax, fourthIndex] = max(relaysSINR(l,:));
                                                        idx = idleRelays(fourthIndex);
                                                        if ~ismember(idx,successfulRelays) && relayState(idx)== 0
                                                            pktStartTXTime(idx) = pktEndTXTime(src) + ST  + aTurnaroudTime + SIFS; % pkt start transmission time
                                                            pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                                            fprintf('\nBase Channel Case: Relay %d of Source %d replies with an Ack',idx,src);
                                                            baseChannel(src) = 0;
                                                            CW(src) = 2;
                                                            NB(src) = 0;
                                                            relayState(idx) = 0;
                                                            successfulRelays(idx) = idx;
                                                            E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                            %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                            relayEnergy(idx) = E1(idx);
                                                            pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                            RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                            break;
                                                        else
                                                            fprintf('\nBase Channel Case: Source %d could not find any relays, "Packet lost"', src);
                                                            baseChannel(src) = 0;
                                                            CW(src) = 2;
                                                            NB(src) = 0;
                                                            pktsLost(src) = pktsLost(src) + 1;
                                                            break;
                                                        end
                                                    end
                                                end
                                            end
                                        else
                                            if length(b) == 5
                                                [val, ind] = max(relaysSINR(l,:));
                                                idx = idleRelays(ind);
                                                if ~ismember(idx,successfulRelays) && relayState(idx)== 0
                                                    pktStartTXTime(idx) = pktEndTXTime(src) + ST  + aTurnaroudTime + SIFS; % pkt start transmission time
                                                    pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                                    fprintf('\nBase Channel Case: Relay %d of Source %d replies with an Ack at %f CW2',idx,src);
                                                    baseChannel(src) = 0;
                                                    CW(src) = 2;
                                                    NB(src) = 0;
                                                    relayState(idx) = 0;
                                                    successfulRelays(idx) = idx;
                                                    E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                    %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                    relayEnergy(idx) = E1(idx);
                                                    pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                    RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                    break;
                                                else
                                                    relaysSINR(l,ind) = -500;
                                                    [secondMax, secondIndex] = max(relaysSINR(l,:));
                                                    idx = idleRelays(secondIndex);
                                                    if ~ismember(idx,successfulRelays) && relayState(idx)== 0
                                                        pktStartTXTime(idx) = pktEndTXTime(src) + ST  + aTurnaroudTime + SIFS; % pkt start transmission time
                                                        pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                                        fprintf('\nBase Channel Case: Relay %d of Source %d replies with an Ack',idx,src);
                                                        baseChannel(src) = 0;
                                                        CW(src) = 2;
                                                        NB(src) = 0;
                                                        relayState(idx) = 0;
                                                        successfulRelays(idx) = idx;
                                                        E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                        %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                        relayEnergy(idx) = E1(idx);
                                                        pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                        RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                        break;
                                                    else
                                                        relaysSINR(l,secondIndex) = -500;
                                                        [thirdMax, thirdIndex] = max(relaysSINR(l,:));
                                                        idx = idleRelays(thirdIndex);
                                                        if ~ismember(idx,successfulRelays) && relayState(idx)== 0
                                                            pktStartTXTime(idx) = pktEndTXTime(src) + ST  + aTurnaroudTime + SIFS; % pkt start transmission time
                                                            pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                                            fprintf('\nBase Channel Case: Relay %d of Source %d replies with an Ack',idx,src);
                                                            baseChannel(src) = 0;
                                                            CW(src) = 2;
                                                            NB(src) = 0;
                                                            relayState(idx) = 0;
                                                            successfulRelays(idx) = idx;
                                                            E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                            %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                            relayEnergy(idx) = E1(idx);
                                                            pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                            RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                            break;
                                                        else
                                                            relaysSINR(l,thirdIndex) = -500;
                                                            [fourthMax, fourthIndex] = max(relaysSINR(l,:));
                                                            idx = idleRelays(fourthIndex);
                                                            if ~ismember(idx,successfulRelays) && relayState(idx)== 0
                                                                pktStartTXTime(idx) = pktEndTXTime(src) + ST  + aTurnaroudTime + SIFS; % pkt start transmission time
                                                                pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                                                fprintf('\nBase Channel Case: Relay %d of Source %d replies with an Ack',idx,src);
                                                                baseChannel(src) = 0;
                                                                CW(src) = 2;
                                                                NB(src) = 0;
                                                                relayState(idx) = 0;
                                                                successfulRelays(idx) = idx;
                                                                E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                                %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                                relayEnergy(idx) = E1(idx);
                                                                pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                                RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                                break;
                                                            else
                                                                fprintf('\nBase Channel Case: Source %d could not find any relays, "Packet lost"', src);
                                                                baseChannel(src) = 0;
                                                                CW(src) = 2;
                                                                NB(src) = 0;
                                                                pktsLost(src) = pktsLost(src) + 1;
                                                                break;
                                                            end
                                                        end
                                                    end
                                                end
                                            else
                                                fprintf('\nBase Channel Case: Relay Number exceeds MaxNumberRelays = 5');
                                                pktsLost(src) = pktsLost(src) + 1;
                                            end%end length(b) == 5
                                        end %end length(b) == 4
                                    end %end length(b) == 3
                                end %end length(b) == 2
                            end %end length(b) == 1
                            break;
                        else
                            csmaTime(src) = csmaTime(src) + TCCA;
                            fprintf('\nBase Channel Case: Node %d performs CCA,%d', src,CW(src));
                            E1(src)=E1(src)+ TCCA*Il*V; % stays listening for TCCA time
                            E3(src)=E3(src)+ TCCA*Il*V; % stays listening for TCCA time
                        end
                    else
                        CW(src) = 2;
                        NB(src) = NB(src) + 1;
                        backoff(i) = min(backoff(i)*BP + 1, macMaxBE*BP);
                        if NB(src) <= macMaxCSMABackoffs
                            backoff(i) = floor(2^backoff(i)*BP*rand);
                            csmaTime(src) = backoff(i)*BP + TCCA;
                            E1(src)=E1(src)+ TCCA*Il*V; % stays listening for TCCA time
                            E3(src)=E3(src)+ TCCA*Il*V; % stays listening for TCCA time
                            fprintf('\nBase Busy Channel Case: NB[%d] = %d', src,NB(src));
                        else
                            fprintf('\nBase Channel Case: Node %d terminates; NB[ %d ] exceeds macMaxCSMABackoffs = 5', src, src);
                            pktsLost(src) = pktsLost(src) + 1;
                            break;
                        end
                    end
                end
                %**********************************************************************************************************************************************************
            else %the source is outside a distance of radius r centered at C
                fprintf('\nReserved Channel Case: Node %d switches onto reserved channel', src);
                CW(src) = 2;
                NB(src) = 0;
                csmaTime(src) = 0;
                pktStartTXTime(src) = 0;
                pktEndTXTime(src) = 0;
                backoff(src) = floor(2^minBE*rand);
                csmaTime(src) = backoff(i)*BP + TCCA;
                E1(src)=E1(src)+ TCCA*Il*V; % stays listening for TCCA time
                E3(src)=E3(src)+ TCCA*Il*V; % stays listening for TCCA time
                fprintf('\nReserved Channel Case: Node %d performs CCA', src);
                while NB(src) <= macMaxCSMABackoffs
                    if baseChannel(src)
                        numberMessagesTransmitted = numberMessagesTransmitted + 1;
                        CW(src) = CW(src) - 1;
                        if CW(src) == 0
                            rc1(src) = 0;
                            CW(src) = 2;
                            %the source sends data packet
                            pktStartTXTime(src) = csmaTime(src) + ST + SIFS; % pkt start transmission time
                            pktEndTXTime(src) = pktStartTXTime(src) + aTurnaroudTime + m*P*TB; % pkt end transmission time
                            fprintf('\nReserved Channel Case: Node %d transmits data packet', src);
                            pktsTransmitted(src) = pktsTransmitted(src) + 1;
                            srcPktsTXed(src) = srcPktsTXed(src) + 1;
                            E1(src)=E1(src)+ (m*P*TB*It*V); %energy cost for transmiting m packets
                            %E1(src)=E1(src)+ (macACKWaitDuration*Il*V); % %energy cost to wait (listening) an ack for the last packet transmitted
                            %E1(src)=E1(src)+ (Ack*TB*Ir*V); % energy cost for receiving one ack
                            E3(src)=E3(src)+ (m*P*TB*It*V); %energy cost for transmiting m packets
                            %E3(src)=E3(src)+ (macACKWaitDuration*Il*V); % %energy cost to wait (listening) an ack for the last packet transmitted
                            %E3(src)=E3(src)+ (Ack*TB*Ir*V);
                            for l = 1:sizeT
                                for k=1:sizeR% size of all relays
                                    s2rdist = topo_dist(transmitters(l),idleRelays(k)); % a source to relays distance    %%% add SINR instead
                                    src2rlyDist(k) = s2rdist;
                                    [pow_received, Path_Loss, snr]  = Implant_PL_Model(s2rdist,N_point,2);
                                    P_Received(l,k) = pow_received;
                                    relayInterferences(l,k) =  sum(P_Received(:,k)) - P_Received(l,k);
                                    relaysSINR(l,k) = P_Received(l,k) - (relayInterferences(l,k) - Noise_Floor);
                                end
                            end
                            %Optimal relay mechanism
                            maxRelaySINR = zeros(1,sizeT);
                            for l = 1:sizeT
                                for k =1:sizeR
                                    relayInterferences(l,k) =  sum(P_Received(:,k)) - P_Received(l,k);
                                    relaysSINR(l,k) = P_Received(l,k) - (relayInterferences(l,k) - Noise_Floor);
                                end
                            end
                            fprintf('\nReserved Channel; Relays SINR;');
                            disp(relaysSINR(1,:));
                            avg = mean(relaysSINR(l,:));
                            s = std(relaysSINR(l,:));
                            b = find(relaysSINR(l,:) >= avg);
                            [val, ind]= max(relaysSINR(l,:));
                            idx = idleRelays(ind);
                            if length(b) == 1
                                [val, ind] = max(relaysSINR(l,:));
                                idx = idleRelays(ind);
                                if ~ismember(idx,reservedCaseRelays) && relayState(idx) == 0 %if idle relay
                                    pktStartTXTime(idx) = pktEndTXTime(src) + ST + 2 * aTurnaroudTime + SIFS; % pkts start transmission time
                                    pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkts end transmission time
                                    fprintf('\nReserved Channel Case: Relay %d of Source %d replies with an Ack',idx,src);
                                    E1(src)=E1(src)+ (m*P*TB*It*V); %energy cost for transmiting m packets
                                    E1(src)=E1(src)+ (macACKWaitDuration*Il*V); % %energy cost to wait (listening) an ack for the last packet transmitted
                                    E1(src)=E1(src)+ (Ack*TB*Ir*V); % energy cost for receiving one ack
                                    CW(src) = 2;
                                    %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                    RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                    reservedCaseRelays(idx) = idx;
                                    relayEnergy(idx) = E1(idx);
                                    NB(src) = 0;
                                    break;
                                else
                                    fprintf('\nReserved Channel Case: Source %d could not find any relays, "Packet lost"', src);
                                    CW(src) = 2;
                                    NB(src) = 0;
                                    pktsLost(src) = pktsLost(src) + 1;
                                    break;
                                end
                            else
                                if length(b) == 2
                                    [val, ind] = max(relaysSINR(l,:));
                                    idx = idleRelays(ind);
                                    if ~ismember(idx,reservedCaseRelays) && relayState(idx) == 0
                                        pktStartTXTime(idx) = pktEndTXTime(src) + ST + 2 * aTurnaroudTime + SIFS; % pkts start transmission time
                                        pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkts end transmission time
                                        fprintf('\nReserved Channel Case: Relay %d of Source %d replies with an Ack',idx,src);
                                        E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                        %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                        relayEnergy(idx) = E1(idx);
                                        %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                        RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                        reservedCaseRelays(idx) = idx;
                                        CW(src) = 2;
                                        NB(src) = 0;
                                        break;
                                    else
                                        relaysSINR(l,ind) = -500;
                                        [secondMax, secondIndex] = max(relaysSINR(l,:));
                                        idx = idleRelays(secondIndex);
                                        if ~ismember(idx,reservedCaseRelays) && relayState(idx) == 0
                                            pktStartTXTime(idx) = pktEndTXTime(src) + ST + 2 * aTurnaroudTime + SIFS; % pkts start transmission time
                                            pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkts end transmission time
                                            fprintf('\nReserved Channel Case: Relay %d of Source %d replies with an Ack',idx,src);
                                            E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                            %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                            relayEnergy(idx) = E1(idx);
                                            %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                            RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                            reservedCaseRelays(idx) = idx;
                                            CW(src) = 2;
                                            NB(src) = 0;
                                            break;
                                        else
                                            fprintf('\nReserved Channel Case: Source %d could not find any relays, "Packet lost"', src);
                                            CW(src) = 2;
                                            NB(src) = 0;
                                            pktsLost(src) = pktsLost(src) + 1;
                                            break;
                                        end
                                    end
                                else
                                    if length(b) == 3
                                        [val, ind] = max(relaysSINR(l,:));
                                        idx = idleRelays(ind);
                                        if ~ismember(idx,reservedCaseRelays) && relayState(idx) == 0
                                            pktStartTXTime(idx) = pktEndTXTime(src) + ST + 2 * aTurnaroudTime + SIFS; % pkts start transmission time
                                            pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkts end transmission time
                                            fprintf('\nReserved Channel Case: Relay %d of Source %d replies with an Ack',idx,src);
                                            E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                            %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                            relayEnergy(idx) = E1(idx);
                                            %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                            RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                            reservedCaseRelays(idx) = idx;
                                            CW(src) = 2;
                                            NB(src) = 0;
                                            break;
                                        else
                                            relaysSINR(l,ind) = -500;
                                            [secondMax, secondIndex] = max(relaysSINR(l,:));
                                            idx = idleRelays(secondIndex);
                                            if ~ismember(idx,reservedCaseRelays) && relayState(idx) == 0
                                                pktStartTXTime(idx) = pktEndTXTime(src) + ST + 2 * aTurnaroudTime + SIFS; % pkts start transmission time
                                                pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkts end transmission time
                                                fprintf('\nReserved Channel Case: Relay %d of Source %d replies with an Ack',idx,src);
                                                reservedCaseRelays(idx) = idx;
                                                E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                relayEnergy(idx) = E1(idx);
                                                %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                CW(src) = 2;
                                                NB(src) = 0;
                                                break;
                                            else
                                                relaysSINR(l,secondIndex) = -500;
                                                [thirdMax, thirdIndex] = max(relaysSINR(l,:));
                                                idx = idleRelays(thirdIndex);
                                                if ~ismember(idx,reservedCaseRelays) && relayState(idx) == 0
                                                    pktStartTXTime(idx) = pktEndTXTime(src) + ST + 2 * aTurnaroudTime + SIFS; % pkts start transmission time
                                                    pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkts end transmission time
                                                    fprintf('\nReserved Channel Case: Relay %d of Source %d replies with an Ack',idx,src);
                                                    reservedCaseRelays(idx) = idx;
                                                    E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                    %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                    relayEnergy(idx) = E1(idx);
                                                    %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                    RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                    CW(src) = 2;
                                                    NB(src) = 0;
                                                    break;
                                                else
                                                    fprintf('\nReserved Channel Case: Source %d could not find any relays, "Packet lost"', src);
                                                    CW(src) = 2;
                                                    NB(src) = 0;
                                                    pktsLost(src) = pktsLost(src) + 1;
                                                    break;
                                                end
                                            end
                                        end
                                    else
                                        if length(b) == 4
                                            [val, ind] = max(relaysSINR(l,:));
                                            idx = idleRelays(ind);
                                            if ~ismember(idx,reservedCaseRelays) && relayState(idx) == 0
                                                pktStartTXTime(idx) = pktEndTXTime(src) + ST  + 2 * aTurnaroudTime + SIFS; % pkt start transmission time
                                                pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                                fprintf('\nReserved Channel Case: Relay %d of Source %d replies with an Ack',idx,src);
                                                reservedCaseRelays(idx) = idx;
                                                E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                relayEnergy(idx) = E1(idx);
                                                %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                CW(src) = 2;
                                                NB(src) = 0;
                                                break;
                                            else
                                                relaysSINR(l,ind) = -500;
                                                [secondMax, secondIndex] = max(relaysSINR(l,:));
                                                idx = idleRelays(secondIndex);
                                                if ~ismember(idx,reservedCaseRelays) && relayState(idx) == 0
                                                    pktStartTXTime(idx) = pktEndTXTime(src) + ST  + 2 * aTurnaroudTime + SIFS; % pkt start transmission time
                                                    pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                                    fprintf('\nReserved Channel Case: Relay %d of Source %d replies with an Ack',idx,src);
                                                    E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                    %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                    relayEnergy(idx) = E1(idx);
                                                    %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                    RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                    reservedCaseRelays(idx) = idx;
                                                    CW(src) = 2;
                                                    NB(src) = 0;
                                                    break;
                                                else
                                                    relaysSINR(l,secondIndex) = -500;
                                                    [thirdMax, thirdIndex] = max(relaysSINR(l,:));
                                                    idx = idleRelays(thirdIndex);
                                                    if ~ismember(idx,reservedCaseRelays) && relayState(idx) == 0
                                                        pktStartTXTime(idx) = pktEndTXTime(src) + ST  + 2 * aTurnaroudTime + SIFS; % pkt start transmission time
                                                        pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                                        fprintf('\nReserved Channel Case: Relay %d of Source %d replies with an Ack',idx,src);
                                                        E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                        %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                        relayEnergy(idx) = E1(idx);
                                                        %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                        RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                        reservedCaseRelays(idx) = idx;
                                                        CW(src) = 2;
                                                        NB(src) = 0;
                                                        break;
                                                    else
                                                        relaysSINR(l,thirdIndex) = -500;
                                                        [fourthMax, fourthIndex] = max(relaysSINR(l,:));
                                                        idx = idleRelays(fourthIndex);
                                                        if ~ismember(idx,reservedCaseRelays) && relayState(idx) == 0
                                                            pktStartTXTime(idx) = pktEndTXTime(src) + ST  + 2 * aTurnaroudTime + SIFS; % pkt start transmission time
                                                            pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                                            fprintf('\nReserved Channel Case: Relay %d of Source %d replies with an Ack ',idx,src);
                                                            E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                            %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                            relayEnergy(idx) = E1(idx);
                                                            %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                            RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                            reservedCaseRelays(idx) = idx;
                                                            CW(src) = 2;
                                                            NB(src) = 0;
                                                            break;
                                                        else
                                                            fprintf('\nReserved Channel Case: Source %d could not find any relays, "Packet lost"', src);
                                                            CW(src) = 2;
                                                            NB(src) = 0;
                                                            pktsLost(src) = pktsLost(src) + 1;
                                                            break;
                                                        end
                                                    end
                                                end
                                            end
                                        else
                                            if length(b) == 5
                                                [val, ind] = max(relaysSINR(l,:));
                                                idx = idleRelays(ind);
                                                if ~ismember(idx,reservedCaseRelays) && relayState(idx) == 0
                                                    pktStartTXTime(idx) = pktEndTXTime(src) + ST  + aTurnaroudTime + SIFS; % pkt start transmission time
                                                    pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                                    fprintf('\nReserved Channel Case: Relay %d of Source %d replies with an Ack',idx,src);
                                                    CW(src) = 2;
                                                    NB(src) = 0;
                                                    E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                    %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                    relayEnergy(idx) = E1(idx);
                                                    reservedCaseRelays(idx) = idx;
                                                    %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                    RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                    break;
                                                else
                                                    relaysSINR(l,ind) = -500;
                                                    [secondMax, secondIndex] = max(relaysSINR(l,:));
                                                    idx = idleRelays(secondIndex);
                                                    if ~ismember(idx,reservedCaseRelays) && relayState(idx) == 0
                                                        pktStartTXTime(idx) = pktEndTXTime(src) + ST  + aTurnaroudTime + SIFS; % pkt start transmission time
                                                        pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                                        fprintf('\nReserved Channel Case: Relay %d of Source %d replies with an Ack',idx,src);
                                                        CW(src) = 2;
                                                        NB(src) = 0;
                                                        reservedCaseRelays(idx) = idx;
                                                        E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                        %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                        relayEnergy(idx) = E1(idx);
                                                        %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                        RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                        break;
                                                    else
                                                        relaysSINR(l,secondIndex) = -500;
                                                        [thirdMax, thirdIndex] = max(relaysSINR(l,:));
                                                        idx = idleRelays(thirdIndex);
                                                        if ~ismember(idx,reservedCaseRelays) && relayState(idx) == 0
                                                            pktStartTXTime(idx) = pktEndTXTime(src) + ST  + aTurnaroudTime + SIFS; % pkt start transmission time
                                                            pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                                            fprintf('\nReserved Channel Case: Relay %d of Source %d replies with an Ack',idx,src);
                                                            CW(src) = 2;
                                                            NB(src) = 0;
                                                            reservedCaseRelays(idx) = idx;
                                                            E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                            %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                            relayEnergy(idx) = E1(idx);
                                                            %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                            RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                            break;
                                                        else
                                                            relaysSINR(l,thirdIndex) = -500;
                                                            [fourthMax, fourthIndex] = max(relaysSINR(l,:));
                                                            idx = idleRelays(fourthIndex);
                                                            if ~ismember(idx,reservedCaseRelays) && relayState(idx) == 0
                                                                pktStartTXTime(idx) = pktEndTXTime(src) + ST  + aTurnaroudTime + SIFS; % pkt start transmission time
                                                                pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                                                fprintf('\nReserved Channel Case: Relay %d of Source %d replies with an Ack',idx,src);
                                                                CW(src) = 2;
                                                                NB(src) = 0;
                                                                reservedCaseRelays(idx) = idx;
                                                                E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                                %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                                relayEnergy(idx) = E1(idx);
                                                                %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                                RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                                break;
                                                            else
                                                                fprintf('\nReserved Channel Case: Source %d could not find any relays, "Packet lost"', src);
                                                                baseChannel(src) = 0;
                                                                CW(src) = 2;
                                                                NB(src) = 0;
                                                                pktsLost(src) = pktsLost(src) + 1;
                                                                break;
                                                            end
                                                        end
                                                    end
                                                end
                                            else
                                                fprintf('\nReserved Channel Case: Relay Number exceeds MaxNumberRelays = 5');
                                                pktsLost(src) = pktsLost(src) + 1;
                                            end%end length (b ) == 5
                                        end%end length(b) == 4
                                    end %end length(b) == 3
                                end%end length(b) == 2
                            end %end length(b) == 1
                            break;
                        else
                            csmaTime(src) = csmaTime(src) + TCCA;
                            E1(src)=E1(src)+ TCCA*Il*V; % stays listening for TCCA time
                            E3(src)=E3(src)+ TCCA*Il*V; % stays listening for TCCA time
                            fprintf('\nReserved Channel Case: Node %d performs CCA,%d',src,CW(src));
                        end
                    else
                        CW(src) = 2;
                        NB(src) = NB(src) + 1;
                        backoff(i) = min(backoff(i) + 1, macMaxBE);
                        if NB(src) <= macMaxCSMABackoffs
                            backoff(i) = floor(2^backoff(i)*BP*rand);
                            csmaTime(src) = backoff(i)*BP + TCCA;
                            E3(src)=E3(src)+ TCCA*Il*V;
                            E1(src)=E1(src)+ TCCA*Il*V; % stays listening for TCCA time
                            fprintf('\nReserved Busy Channel Case: NB[%d] = %d', src,NB(src));
                        else
                            fprintf('\nReserved Channel Case: Node %d terminates; NB[ %d ] exceeds macMaxCSMABackoffs = 5', src, src);
                            pktsLost(src) = pktsLost(src) + 1;
                            break;
                        end
                    end
                end % end while NB
            end %end d <  r
            %*************************************************************************
            %********************************
        else %not member (i.e  interfering source)
            fprintf('\nInterfering Case CW2: Node %d is an interfering source, extends the contention window', src);
            u = sortedTransmitters(sizeSB);
            pktEndTXTime(u) = maxBO * BP + 2 * TCCA + SIFS + 2* ST  + aTurnaroudTime + m * P * TB + SIFS + Ack * TB + SIFS + macACKWaitDuration;
            previousfinishingTime = maxBO * BP + TCCA + SIFS + ST + aTurnaroudTime + m * P * TB + SIFS + Ack * TB + SIFS + macACKWaitDuration;
            csmaTime(src) = previousfinishingTime + SIFS;
            csmaTime(src) = csmaTime(src) + backoff(i) * BP + TCCA;
            fprintf('\nInterfering Case: Node %d performs CCA CW2', src);
            NB(src) = 0;
            CW(src) = 2;
            E3(src)=E3(src)+ TCCA*Il*V; % stays listening for TCCA time
            E1(src)=E1(src)+ TCCA*Il*V; % stays listening for TCCA time
            %determine the channel based on new backoff
            [c,d,indSrcinIS] = intersect(src,sortedInterferingSources);
            srcBKvalue = ISBackoff(indSrcinIS);
            for cc = 1 : sizeIS
                if srcBKvalue == ISBackoff(cc) && cc ~= indSrcinIS
                    ISBaseChannel(indSrcinIS) = 0; % busy channel
                else
                    ISBaseChannel(indSrcinIS) = 1; % clear channel
                end
            end
            if ISBaseChannel(indSrcinIS) == 0 % busy channel
                fprintf('\nInterfering Channel Case: %d "Busy Channel"',ISBaseChannel(indSrcinIS))
            else
                fprintf('\nInterfering Channel Case: %d "Clear Channel"',ISBaseChannel(indSrcinIS))
            end
            while NB(src) <= macMaxCSMABackoffs
                if ISBaseChannel(indSrcinIS) == 1 %clear channel
                    CW(src) = CW(src) - 1;
                    if CW(src) == 0
                        numberMessagesTransmitted = numberMessagesTransmitted + 1;
                        CW(src) = 2;
                        NB(src) = 0;
                        baseChannel(src) = 0;
                        ISBaseChannel(indSrcinIS) = 0;
                        pktStartTXTime(src) = csmaTime(src) + ST + SIFS; % pkt start transmission time
                        pktEndTXTime(src) = pktStartTXTime(src) + aTurnaroudTime + m * P * TB; % pkt end transmission time
                        fprintf('\nInterfering Case: Node %d transmits data packet at %f CW2', src,pktStartTXTime(src));
                        E1(src)=E1(src)+ (m*P*TB*It*V); %energy cost for transmiting m packets
                        %E1(src)=E1(src)+ (macACKWaitDuration*Il*V); % %energy cost to wait (listening) an ack for the last packet transmitted
                        %E1(src)=E1(src)+ (Ack*TB*Ir*V); % energy cost for receiving one ack
                        %E3(src)=E3(src)+ (m*P*TB*It*V); %energy cost for transmiting m packets
                        %E3(src)=E3(src)+ (macACKWaitDuration*Il*V); % %energy cost to wait (listening) an ack for the last packet transmitted
                        %E3(src)=E3(src)+ (Ack*TB*Ir*V); % energy cost for receiving one ack
                        pktsTransmitted(src) = pktsTransmitted(src) + 1;
                        srcPktsTXed(src) = srcPktsTXed(src) + 1;
                        for l = 1:sizeT
                            for k=1:sizeR% size of all relays
                                s2rdist = topo_dist(transmitters(l),idleRelays(k)); % a source to relays distance    %%% add SINR instead
                                src2rlyDist(k) = s2rdist;
                                [pow_received, Path_Loss, snr]  = Implant_PL_Model(s2rdist,N_point,2);
                                P_Received(l,k) = pow_received;
                                relayInterferences(l,k) =  sum(P_Received(:,k)) - P_Received(l,k);
                                relaysSINR(l,k) = P_Received(l,k) - Noise_Floor; %- (relayInterferences(l,k)
                            end
                        end
                        %Optimal relay mechanism
                        maxRelaySINR = zeros(1,sizeT);
                        for l = 1:sizeT
                            for k =1:sizeR
                                relayInterferences(l,k) =  sum(P_Received(:,k)) - P_Received(l,k);
                                relaysSINR(l,k) = P_Received(l,k) - Noise_Floor; %- (relayInterferences(l,k)%
                            end
                        end
                        fprintf('\nInterfering Case; Relays SINR;');
                        disp(relaysSINR(1,:));
                        avg = mean(relaysSINR(l,:));
                        s = std(relaysSINR(l,:));
                        b = find(relaysSINR(l,:) >= sinrThr);
                        [val, ind]= max(relaysSINR(l,:));
                        idx = idleRelays(ind);
                        value = ismember(idx,successfulRelays);
                        if length(b) == 1
                            [val, ind]= max(relaysSINR(l,:));
                            idx = idleRelays(ind);
                            if ~ismember(idx,interCaseRelays) && relayState(idx) == 0 %if idle relay
                                pktStartTXTime(idx) = pktEndTXTime(src) + ST  + 2 * aTurnaroudTime + SIFS; % pkt start transmission time
                                pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                fprintf('\nInterfering Case CW2: Relay %d of Source %d replies with an Ack at %f CW2',idx,src,pktStartTXTime(idx));
                                interCaseRelays(idx) = idx;
                                relayState(idx) = 0;
                                baseChannel(src) = 0;
                                ISBaseChannel(indSrcinIS) = 0;
                                CW(src) = 2;
                                NB(src) = 0;
                                E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                relayEnergy(idx) = E1(idx);
                                %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                break;
                            else
                                fprintf('\nInterfering Case CW2: Source %d could not find any relays, "Packet lost"', src);
                                baseChannel(src) = 0;
                                ISBaseChannel(indSrcinIS) = 0;
                                CW(src) = 2;
                                NB(src) = 0;
                                pktsLost(src) = pktsLost(src) + 1;
                                break;
                            end
                        else
                            if length(b) == 2
                                [val, ind]= max(relaysSINR(l,:));
                                idx = idleRelays(ind);
                                if ~ismember(idx,interCaseRelays) && relayState(idx) == 0 %if idle relay
                                    pktStartTXTime(idx) = pktEndTXTime(src) + ST  + 2 * aTurnaroudTime + SIFS; % pkt start transmission time
                                    pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                    fprintf('\nInterfering Case CW2: Relay %d of Source %d replies with an Ack at %f CW2',idx,src,pktStartTXTime(idx));
                                    interCaseRelays(idx) = idx;
                                    relayState(idx) = 0;
                                    baseChannel(src) = 0;
                                    relayEnergy(idx) = E1(idx);
                                    ISBaseChannel(indSrcinIS) = 0;
                                    CW(src) = 2;
                                    NB(src) = 0;
                                    E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                    %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                    %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                    RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                    break;
                                else
                                    relaysSINR(l,ind) = -500;
                                    [secondMax, secondIndex] = max(relaysSINR(l,:));
                                    idx = idleRelays(secondIndex);
                                    if ~ismember(idx ,interCaseRelays) && relayState(idx)== 0
                                        pktStartTXTime(idx) = pktEndTXTime(src) + ST  + 2 * aTurnaroudTime + SIFS; % pkt start transmission time
                                        pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                        fprintf('\nInterfering Case CW2: Relay %d of Source %d replies with an Ack at %f CW2',idx,src,pktStartTXTime(idx));
                                        interCaseRelays(idx) = idx;
                                        relayState(idx) = 0;
                                        baseChannel(src) = 0;
                                        ISBaseChannel(indSrcinIS) = 0;
                                        relayEnergy(idx) = E1(idx);
                                        CW(src) = 2;
                                        NB(src) = 0;
                                        E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                        %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                        %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                        RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                        break;
                                    else
                                        fprintf('\nInterfering Case CW2: Source %d could not find any relays, "Packet lost"', src);
                                        baseChannel(src) = 0;
                                        ISBaseChannel(indSrcinIS) = 0;
                                        CW(src) = 2;
                                        NB(src) = 0;
                                        pktsLost(src) = pktsLost(src) + 1;
                                        break;
                                    end
                                end
                            else
                                if length(b) == 3
                                    [val, ind]= max(relaysSINR(l,:));
                                    idx = idleRelays(ind);
                                    if ~ismember(idx,interCaseRelays) && relayState(idx)== 0 %if idle relay
                                        pktStartTXTime(idx) = pktEndTXTime(src) + ST  + 2 * aTurnaroudTime + SIFS; % pkt start transmission time
                                        pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                        fprintf('\nInterfering Case CW2: Relay %d of Source %d replies with an Ack at %f CW2',idx,src,pktStartTXTime(idx));
                                        interCaseRelays(idx) = idx;
                                        relayState(idx) = 0;
                                        baseChannel(src) = 0;
                                        ISBaseChannel(indSrcinIS) = 0;
                                        relayEnergy(idx) = E1(idx);
                                        CW(src) = 2;
                                        NB(src) = 0;
                                        E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                        %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                        %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                        RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                        break;
                                    else
                                        relaysSINR(l,ind) = -500;
                                        [secondMax, secondIndex] = max(relaysSINR(l,:));
                                        idx = idleRelays(secondIndex);
                                        if ~ismember(idx,interCaseRelays) && relayState(idx)== 0
                                            idx = idleRelays(secondIndex);
                                            pktStartTXTime(idx) = pktEndTXTime(src) + ST  + 2 * aTurnaroudTime + SIFS; % pkt start transmission time
                                            pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                            fprintf('\nInterfering Case CW2: Relay %d of Source %d replies with an Ack at %f CW2',idx,src,pktStartTXTime(idx));
                                            interCaseRelays(idx) = idx;
                                            relayState(idx) = 0;
                                            baseChannel(src) = 0;
                                            ISBaseChannel(indSrcinIS) = 0;
                                            relayEnergy(idx) = E1(idx);
                                            CW(src) = 2;
                                            NB(src) = 0;
                                            E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                            %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                            %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                            RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                            break;
                                        else
                                            relaysSINR(l,secondIndex) = -500;
                                            [thirdMax, thirdIndex] = max(relaysSINR(l,:));
                                            idx = idleRelays(thirdIndex);
                                            value = ismember(idx,successfulRelays);
                                            if ~ismember(idx,interCaseRelays) && relayState(idx)== 0
                                                pktStartTXTime(idx) = pktEndTXTime(src) + ST  + 2 * aTurnaroudTime + SIFS; % pkt start transmission time
                                                pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                                fprintf('\nInterfering Case CW2: Relay %d of Source %d replies with an Ack at %f CW2',idx,src,pktStartTXTime(idx));
                                                interCaseRelays(idx) = idx;
                                                relayState(idx) = 0;
                                                baseChannel(src) = 0;
                                                ISBaseChannel(indSrcinIS) = 0;
                                                CW(src) = 2;
                                                NB(src) = 0;
                                                relayEnergy(idx) = E1(idx);
                                                E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                break;
                                            else
                                                fprintf('\nInterfering Case CW2: Source %d could not find any relays, "Packet lost"', src);
                                                baseChannel(src) = 0;
                                                ISBaseChannel(indSrcinIS) = 0;
                                                CW(src) = 2;
                                                NB(src) = 0;
                                                pktsLost(src) = pktsLost(src) + 1;
                                                break;
                                            end
                                        end
                                    end
                                else
                                    if length(b) == 4
                                        [val, ind]= max(relaysSINR(l,:));
                                        idx = idleRelays(ind);
                                        if ~ismember(idx,interCaseRelays) && relayState(idx)== 0
                                            pktStartTXTime(idx) = pktEndTXTime(src) + ST  + 2 * aTurnaroudTime + SIFS; % pkt start transmission time
                                            pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                            fprintf('\nInterfering Case CW2: Relay %d of Source %d replies with an Ack at %f CW2',idx,src,pktStartTXTime(idx));
                                            interCaseRelays(idx) = idx;
                                            relayState(idx) = 0;
                                            relayState(idx) = 0;
                                            relayEnergy(idx) = E1(idx);
                                            baseChannel(src) = 0;
                                            ISBaseChannel(indSrcinIS) = 0;
                                            CW(src) = 2;
                                            NB(src) = 0;
                                            E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                            %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                            %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                            RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                            break;
                                        else
                                            relaysSINR(l,ind) = -500;
                                            [secondMax, secondIndex] = max(relaysSINR(l,:));
                                            idx = idleRelays(secondIndex);
                                            if ~ismember(idx,interCaseRelays) && relayState(idx)== 0
                                                pktStartTXTime(idx) = pktEndTXTime(src) + ST  + 2 * aTurnaroudTime + SIFS; % pkt start transmission time
                                                pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                                fprintf('\nInterfering Case CW2: Relay %d of Source %d replies with an Ack at %f CW2',idx,src,pktStartTXTime(idx));
                                                interCaseRelays(idx) = idx;
                                                relayState(idx) = 0;
                                                baseChannel(src) = 0;
                                                ISBaseChannel(indSrcinIS) = 0;
                                                CW(src) = 2;
                                                NB(src) = 0;
                                                relayEnergy(idx) = E1(idx);
                                                E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                % E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                break;
                                            else
                                                relaysSINR(l,secondIndex) = -500;
                                                [thirdMax, thirdIndex] = max(relaysSINR(l,:));
                                                idx = idleRelays(thirdIndex);
                                                if ~ismember(idx,interCaseRelays) && relayState(idx)== 0
                                                    pktStartTXTime(idx) = pktEndTXTime(src) + ST  + 2 * aTurnaroudTime + SIFS; % pkt start transmission time
                                                    pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                                    fprintf('\nInterfering Case CW2: Relay %d of Source %d replies with an Ack at %f CW2',idx,src,pktStartTXTime(idx));
                                                    interCaseRelays(idx) = idx;
                                                    relayState(idx) = 0;
                                                    baseChannel(src) = 0;
                                                    relayEnergy(idx) = E1(idx);
                                                    ISBaseChannel(indSrcinIS) = 0;
                                                    CW(src) = 2;
                                                    NB(src) = 0;
                                                    E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                    %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                    %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                    RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                    break;
                                                else
                                                    relaysSINR(l,thirdIndex) = -500;
                                                    [fourthMax, fourthIndex] = max(relaysSINR(l,:));
                                                    idx = idleRelays(fourthIndex);
                                                    if ~ismember(idx,interCaseRelays) && relayState(idx)== 0
                                                        pktStartTXTime(idx) = pktEndTXTime(src) + ST  + 2* aTurnaroudTime + SIFS; % pkt start transmission time
                                                        pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                                        fprintf('\nInterfering Case CW2: Relay %d of Source %d replies with an Ack at %f CW2',idx,src,pktStartTXTime(idx));
                                                        interCaseRelays(idx) = idx;
                                                        relayState(idx) = 0;
                                                        baseChannel(src) = 0;
                                                        ISBaseChannel(indSrcinIS) = 0;
                                                        relayEnergy(idx) = E1(idx);
                                                        CW(src) = 2;
                                                        NB(src) = 0;
                                                        E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                        %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                        %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                        RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                        break;
                                                    else
                                                        fprintf('\nInterfering Case CW2: Source %d could not find any relays, "Packet lost"', src);
                                                        baseChannel(src) = 0;
                                                        ISBaseChannel(indSrcinIS) = 0;
                                                        CW(src) = 2;
                                                        NB(src) = 0;
                                                        pktsLost(src) = pktsLost(src) + 1;
                                                        break;
                                                    end
                                                end
                                            end
                                        end
                                    else
                                        if length(b) == 5
                                            [val, ind]= max(relaysSINR(l,:));
                                            idx = idleRelays(ind);
                                            if ~ismember(ind,interCaseRelays) && relayState(idx)== 0 %if idle relay
                                                pktStartTXTime(idx) = pktEndTXTime(src) + ST  + aTurnaroudTime + SIFS; % pkt start transmission time
                                                pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                                fprintf('\nInterfering Case: Relay %d of Source %d replies with an Ack at %f CW2',idx,src,pktStartTXTime(idx));
                                                interCaseRelays(idx) = idx;
                                                relayState(idx) = 0;
                                                baseChannel(src) = 0;
                                                relayEnergy(idx) = E1(idx);
                                                CW(src) = 2;
                                                NB(src) = 0;
                                                E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                % E1(idx)=E1(idx) + Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                break;
                                            else
                                                relaysSINR(l,ind) = -500;
                                                [secondMax, secondIndex] = max(relaysSINR(l,:));
                                                idx = idleRelays(secondIndex);
                                                if ~ismember(idx,interCaseRelays) && relayState(idx)== 0
                                                    pktStartTXTime(idx) = pktEndTXTime(src)+ ST + aTurnaroudTime + SIFS; % pkt start transmission time
                                                    pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                                    fprintf('\nInterfering Case: Relay %d of Source %d replies with an Ack at %f CW2',idx,src,pktStartTXTime(idx));
                                                    interCaseRelays(idx) = idx;
                                                    relayState(idx) = 0;
                                                    baseChannel(src) = 0;
                                                    relayEnergy(idx) = E1(idx);
                                                    CW(src) = 2;
                                                    NB(src) = 0;
                                                    value = ismember(idx,successfulRelays);
                                                    E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                    %E1(idx)=E1(idx) + Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                    %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                    RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                    break;
                                                else
                                                    relaysSINR(l,secondIndex) = -500;
                                                    [thirdMax, thirdIndex] = max(relaysSINR(l,:));
                                                    idx = idleRelays(thirdIndex);
                                                    if ~ismember(idx,interCaseRelays) && relayState(idx)== 0
                                                        pktStartTXTime(idx) = pktEndTXTime(src) + ST  + aTurnaroudTime + SIFS; % pkt start transmission time
                                                        pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                                        fprintf('\nInterfering Case: Relay %d of Source %d replies with an Ack at %f CW2',idx,src,pktStartTXTime(idx));
                                                        interCaseRelays(idx) = idx;
                                                        relayState(idx) = 0;
                                                        baseChannel(src) = 0;
                                                        relayEnergy(idx) = E1(idx);
                                                        CW(src) = 2;
                                                        NB(src) = 0;
                                                        E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                        %E1(idx)=E1(idx)+Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                        %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                        RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                        break;
                                                    else
                                                        relaysSINR(l,thirdIndex) = -500;
                                                        [fourthMax, fourthIndex] = max(relaysSINR(l,:));
                                                        idx = idleRelays(fourthIndex);
                                                        if ~ismember(idx,interCaseRelays) && relayState(idx)== 0
                                                            pktStartTXTime(idx) = pktEndTXTime(src) + ST  + aTurnaroudTime + SIFS; % pkt start transmission time
                                                            pktEndTXTime(idx)= pktStartTXTime(idx) + Ack*TB + SIFS; % pkt end transmission time
                                                            fprintf('\nInterfering Case: Relay %d of Source %d replies with an Ack at %f CW2',idx,src,pktStartTXTime(idx));
                                                            interCaseRelays(idx) = idx;
                                                            relayState(idx) = 0;
                                                            baseChannel(src) = 0;
                                                            relayEnergy(idx) = E1(idx);
                                                            CW(src) = 2;
                                                            NB(src) = 0;
                                                            E1(idx)=E1(idx)+ m*P*TB*Ir*V; % energy for receiving m pkts
                                                            %E1(idx)=E1(idx)+ Ack*TB*It*V; % energy cost for transmitting an ack to source
                                                            %pktsTransmitted(idx) = pktsTransmitted(idx) + 1;
                                                            RpktsTransmitted(idx) = RpktsTransmitted(idx) + 1;
                                                            break;
                                                        else
                                                            fprintf('\nInterfering Case: Source %d could not find any relays, "Packet lost"', src);
                                                            baseChannel(src) = 0;
                                                            CW(src) = 2;
                                                            NB(src) = 0;
                                                            pktsLost(src) = pktsLost(src) + 1;
                                                            break;
                                                        end
                                                    end
                                                end
                                            end
                                        else
                                            fprintf('\nInterfering Channel Case: Relay Number exceeds MaxNumberRelays = 5');
                                            pktsLost(src) = pktsLost(src) + 1;
                                        end%end length (b ) == 5
                                    end%end length (b ) == 4
                                end %end length (b ) == 3
                            end %end length (b ) == 2
                        end%end length (b ) == 1
                        break;
                    else
                        csmaTime(src) = csmaTime(src) + TCCA;
                        E1(src)=E1(src)+ TCCA*Il*V; % stays listening for TCCA time
                        E3(src)=E3(src)+ TCCA*Il*V; % stays listening for TCCA time
                        fprintf('\nInterfering Case: Node %d performs CCA CW2,%d',src,CW(src));
                    end
                else
                    CW(src) = 2;
                    NB(src) = NB(src) + 1;
                    backoff(i) = min(backoff(i)*BP + 1, macMaxBE*BP);
                    if NB(src) <= macMaxCSMABackoffs
                        backoff(i) = floor(2^backoff(i)*BP*rand);
                        csmaTime(src) = backoff(i)*BP + TCCA;
                        E1(src)=E1(src)+ TCCA*Il*V; % stays listening for TCCA time
                        E3(src)=E3(src)+ TCCA*Il*V; % stays listening for TCCA time
                        fprintf('\nInterfering Channel Case: NB[%d] = %d', src,NB(src));
                    else
                        fprintf('\nInterfering Case: Node %d terminates; NB[ %d ] exceeds macMaxCSMABackoffs = 5', src, src);
                        pktsLost(src) = pktsLost(src) + 1;
                        break;
                    end
                end
            end %endwhile(NB)
        end %end ismember
    end%end for

    counter = counter + 1;

    %TDMA starts
    succRelays = find(successfulRelays > 0);
    sizeSR = length(succRelays);
    resRelays = find(reservedCaseRelays > 0);
    sizeRR = length(resRelays);
    interRelays = find(interCaseRelays > 0);
    sizeIR = length(interRelays);
    sizeSlots =  sizeSR + sizeRR  + sizeIR;
    fprintf('\n/***TDMA SCHEDULE START OF %d SLOTS***/\n',sizeSlots)

    PETXTsuccRelays = pktEndTXTime(succRelays);
    PETXTresRelays = pktEndTXTime(resRelays);
    PETXTinterRelays = pktEndTXTime(interRelays);

    sortedPETXTsucc = sort(PETXTsuccRelays);

    sortedPETXTres = sort(PETXTresRelays);
    sortedPETXTinter = sort(PETXTinterRelays);
    slotsAssigned = zeros(1,sizeSlots);
    concatenatedrelays = [succRelays, resRelays, interRelays];
    sortedconcatenatedrelays = sort(concatenatedrelays);
    [C, ia, ic] = unique(sortedconcatenatedrelays);
    sizeRelays = length(C);
    relaySlotAssigned = zeros(1,sizeRelays);

    for t = 1:sizeRelays
        for z = 1:sizeSlots
            if C(t) == sortedconcatenatedrelays(z)
                relaySlotAssigned(t) = relaySlotAssigned(t) + 1;
            end
        end
    end
    for k = 1 : sizeRelays
        fprintf('Relay %d is assigned %d slots\n',C(k), relaySlotAssigned(k));
        indx = C(k);
        %pktsTransmitted(indx) =  pktsTransmitted(indx) + relaySlotAssigned(k);
        E1(indx)=E1(indx)+ relaySlotAssigned(k) * (m*P*TB*It*V); % energy for transmitting m pkts
        %E1(indx)=E1(indx)+ relaySlotAssigned(k) * (macACKWaitDuration*Il*V + SE);
        %E1(indx)=E1(indx)+ relaySlotAssigned(k) * (Ack*TB*Ir*V); % energy cost for receiving an ack from coordinator
    end
    fprintf('\n/***TDMA SCHEDULE END OF %d SLOTS***/\n',sizeSlots);
    if counter == 1000
        break;
    end
    %     if numberMessagesTransmitted ==100
    %         break;
    %     end
end%end while
srcPktsTXed(4) = 1;
srcPktsTXed(5) = 1;
srcPktsTXed(9) = 1;
srcPktsTXed(10) = 1;
PLR = zeros(1,N);
for n = 2 : N
    PLR(n) = pktsLost(n) / ( pktsLost(n) + srcPktsTXed(n))*100;
end


%ZIGBEE MAC
%**************************************************************************
%node = [0.5 1.2; 0.4 1.9; 0.2 0.5; 0.4 0.6; 0.3 1.4; 0.2 1.9; 0.6 0.3; 0.7 0.8; 0.55 0.9; 0.5 1.6; 0.6 1.9]; %10 nodes
% plot(node(:,1),node(:,2), 'b*')
% hold
% grid on
for j=1:N
    D(1,j) =  topo_dist(1, j);
end
srcPktsT = zeros(1,N);
E2 = zeros(1,N);
counter = 1;
macMaxCSMABackoffs = 4;
numberMessagesTr = 0;
pktsTr = zeros(1,N);
pktsLo = zeros(1,N);
while counter < 1000
    Isources = [];
    NB = zeros(1,N);
    csmaTime = zeros(1,N);
    pktStartTXTime = zeros(1,N);
    pktEndTXTime = zeros(1,N);
    CW = [2 2 2 2 2 2 2 2 2 2];
    sourceState = floor(3*rand(1,N));
%     sourceState(1) = 20;
%     sourceState(2) = 1;
%     sourceState(3) = 1;
%     sourceState(4) = 0;
%     sourceState(5) = 0;%relay
%     sourceState(6) = 1;%source
%     sourceState(7) = 1;
%     sourceState(8) = 1;
%     sourceState(9) = 0;
%     sourceState(10) = 0;
%     sourceState(11) = 1;
%     sourceState(1) = 20;
    sources = find(sourceState>=1 & sourceState < 20);
    S = length(sources); %sources size
    baseChannel = zeros(1,S);
    backoff = floor(2^minBE*rand(1,S));
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
    fprintf('\nZigbee Case: Sources = ');
    disp(sources);
    fprintf('\nZigbee Case: Backoff = ');
    disp(backoff);
    src2Cdist = sortMatrix(:,3)';
    %**************************************************************************
    %number of transmitters on basechannel;
    trans = find(tr>0);
    sizeT = length(trans);
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

    fprintf('Zigbee Case: Sorted backoff: ');
    sortedBackoff = backoff(trans);
    disp(sortedBackoff);
    D = D * 100;
    sizeSB = length(sortedBackoff);
    %for all transmitters set the channel as clear

    intSrcs = setxor(sources,transmitters); %complement of transmitters in sources

    baseChannel(sortedTransmitters)= 1;
    baseChannel(intSrcs ) = 0;

    fprintf('\nZigbee Case channel : ');
    disp(baseChannel);
    minBE = 4;
    %initialize the sources backoff array
    for g =1:S
        csmaTime(g) = backoff(g);
    end
    %Initialize sensors energy
    for i = 1:N
        if r == 1 % coordinator
            E2(i)= 0;
        else
            E2(i)=E2(i)+(2*SE);
        end
    end
    for i = 1:S % for all sources
        src = sources(i);
        NB(src) = 0;
        CW(src) = 2;
        csmaTime(src) = backoff(i)*BP + TCCA;
        E2(src)=E2(src)+ TCCA*Il*V; % stays listening for TCCA time
        E2(src)=E2(src)+ SE;
        fprintf('\nZigbee Case: Node %d performs CCA', src);
        while NB(src) <= macMaxCSMABackoffs
            if baseChannel(src) == 1
                CW(src) = CW(src) - 1;
                if CW(src) == 0
                    %Source transmits m data packets
                    pktStartTXTime(src) = csmaTime(src) + ST + SIFS; % pkts start transmission time
                    pktEndTXTime(src) = pktStartTXTime(src) + m*P*TB; % pkts end transmission time
                    fprintf('\nZigbee Case: Node %d transmits data packet', src);
                    E2(src)=E2(src) + 2*m*SIFS*Il*V + m*P*TB*It*V + m*Ack*TB*Ir*V + 2*(m+1)*SE + m*(macACKWaitDuration*Il*V);
                    pktsTr(src) = pktsTr(src) + 1;
                    numberMessagesTr = numberMessagesTr + 1;
                    NB(src) = 0;
                    CW(src) = 2;
                    srcPktsT(src) = srcPktsT(src) + 1;
                    break;
                else
                    csmaTime(src) = csmaTime(src) + TCCA;
                    fprintf('\nZigbee Case: Node %d performs CCA,%d', src,CW(src));
                    E2(src)=E2(src)+ TCCA*Il*V; % stays listening for TCCA time
                end
            else
                CW(src) = 2;
                NB(src) = NB(src) + 1;
                backoff(i) = min(backoff(i)*BP + 1, macMaxBE*BP);
                if NB(src) <= macMaxCSMABackoffs
                    backoff(i) = floor(2^backoff(i)*BP*rand);
                    csmaTime(src) = backoff(i)*BP + TCCA;
                    E2(src)=E2(src)+ TCCA*Il*V; % stays listening for TCCA time
                    fprintf('\nZigbee Busy Channel Case: NB[%d] = %d', src,NB(src));
                else
                    fprintf('\nZigbee Case: Node %d terminates; NB[ %d ] exceeds macMaxCSMABackoffs = 5', src, src);
                    pktsLo(src) = pktsLo(src) + 1;
                    break;
                end
            end
        end
    end%end for
    counter = counter + 1;
    if counter == 1000
        break;
    end
    if numberMessagesTr == 1000
        break;
    end
end%end while
plr = zeros(1,N);
for n = 2 : N
    plr(n) = pktsLo(n) / ( pktsLo(n) + srcPktsT(n))*100;
end


%Plot Results
%Energy Consumption plot for WBAN versus Zigbee MAC
% plot(E1,'r:+')  % E3 energy of IAA sources only
% hold
% plot( E2,'b:*')
% legend('Energy Consumption "IAA Algo"','Energy Consumption "ZigbeeMAC"')
% xlabel('Node Index');
% ylabel('Energy Consumption [milliJoule]');
% grid on

%Plot of Number of packets lost with IAA algorithm versus ZIGBE
% plot(pktsLost,'r:+')
% hold
% plot(pktsLo, 'b:*')
% legend('IAApktsLost','ZIGBEEpktsLost')
% xlabel('Node Index');
% ylabel('Number of Messages');
% grid on
%
plot(PLR,'r:+')
hold
plot(plr, 'b:*')
legend('Packet Loss Rate "IAA"','Packet Loss Rate "ZIGBEEMAC"')
xlabel('Node Index');
ylabel('Number of Messages');
grid on

% plot(pktsTransmitted,'r:+')
% hold
% plot( pktsLost, 'b:*')
% legend('IAApktsTransmitted','IAApktsLost')
% xlabel('Node Index');
% ylabel('Number of Messages');
% grid on

%**************************************************************************


