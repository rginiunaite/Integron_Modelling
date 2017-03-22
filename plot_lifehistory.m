clear
clf
data = [];

time=[0,1,2];
genFirstPos = [1,3,0];
integrase = [1,0,0];
replication = [0,1,0];
life = [1,1,0];

data(:,1) = time;
data(:,2) = genFirstPos;
data(:,3) = integrase;
data(:,4) = replication;
data(:,5) = life;
S = 100; % size of a circle



for i = 1:(length(time))
  
 
    %check if the cell did not die
    if life(i) == 0
        scatter(time(i),1,'x');
        break
    end
    
    % colour the circle based on which gene is in the first position
    hold on

    % Choose colour according to which gene is in first position
    if genFirstPos(i) == 1
     lineColour = 'blue';
    elseif genFirstPos(i) == 2
     lineColour = 'red';
    elseif genFirstPos(i) == 3
     lineColour = [0.9290    0.6940    0.1250];
    else
    lineColour = [0.4940    0.1840    0.5560];
    end
    
    % Choose line type according to if the integrase is functional or not
    if integrase(i) == 1
        lineType = '-';
    else
        lineType = '--';
    end
    
    plot([time(i), time(i+1)],[1,1],'color',lineColour,'lineStyle',lineType);
    
    % Put a dot if the cell replicated
    if replication(i) == 1
        scatter(time(i),1,S,'MarkerFaceColor',lineColour);
    end
end
shg

