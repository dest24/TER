%ZIGBEE MAC
%**************************************************************************
%node = [0.5 1.2; 0.4 1.9; 0.2 0.5; 0.4 0.6; 0.3 1.4; 0.2 1.9; 0.6 0.3; 0.7 0.8; 0.55 0.9; 0.5 1.6; 0.6 1.9]; %10 nodes
node = [0.5 1.2; rand 2*rand; rand 2*rand; rand 2*rand; rand 2*rand; rand 2*rand; rand 2*rand; rand 2*rand; rand 2*rand; rand 2*rand; rand 2*rand];
% plot(node(:,1),node(:,2), 'b*')
% hold
% grid on
N = 11;
E2 = zeros(1,N);

for j=1:N
    D(1,j) =  topo_dist(1, j);
end
srcPktsT = zeros(1,N);
counter = 1;
macMaxCSMABackoffs = 4;
numberMessagesTr = 0;
pktsTr = zeros(1,N);
pktsLo = zeros(1,N);
for q = 1:20
    while counter < 10
        Isources = [];
        NB = zeros(1,N);
        csmaTime = zeros(1,N);
        pktStartTXTime = zeros(1,N);
        pktEndTXTime = zeros(1,N);
        CW = [2 2 2 2 2 2 2 2 2 2];
        sourceState = floor(5*rand(1,N));
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

        for i = 1:S % for all sources
            src = sources(i);
            NB(src) = 0;
            CW(src) = 2;
            csmaTime(src) = backoff(i)*BP + TCCA;
            fprintf('\nZigbee Case: Node %d performs CCA', src);
            while NB(src) <= macMaxCSMABackoffs
                %if baseChannel(src) == 1
                w = rand();
                if (w > (q / 100))
                    CW(src) = CW(src) - 1;
                    if CW(src) == 0
                        %Source transmits m data packets
                        pktStartTXTime(src) = csmaTime(src) + ST + SIFS; % pkts start transmission time
                        pktEndTXTime(src) = pktStartTXTime(src) + m*P*TB; % pkts end transmission time
                        fprintf('\nZigbee Case: Node %d transmits data packet', src);
                        E2(src)=E2(src)+ m*P*TB*It*V + m*Ack*TB*Ir*V + 2*(m+1)*SE + m*(macACKWaitDuration*Il*V);
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
    end%end while

end
E2(1) = 0;
plot( cumsum(E2),'b:*')
legend('Energy Consumption "IAA Algo"','Energy Consumption "ZigbeeMAC"')
xlabel('Node Index');
ylabel('Energy Consumption [milliJoule]');
grid on




Plot Results
Energy Consumption plot for WBAN versus Zigbee MAC


Plot of Number of packets lost with IAA algorithm versus ZIGBE
plot(pktsLost,'r:+')
hold
plot(pktsLo, 'b:*')
legend('IAApktsLost','ZIGBEEpktsLost')
xlabel('Node Index');
ylabel('Number of Messages');
grid on

plot(PLR,'r:+')
hold
plot(plr, 'b:*')
legend('Packet Loss Rate "IAA"','Packet Loss Rate "ZIGBEEMAC"')
xlabel('Node Index');
ylabel('Number of Messages');
grid on

plot(pktsTransmitted,'r:+')
hold
plot( pktsLost, 'b:*')
legend('IAApktsTransmitted','IAApktsLost')
xlabel('Node Index');
ylabel('Number of Messages');
grid on

%**************************************************************************


