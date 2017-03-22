clear

data = [];

time=[1,2,3];
genFirstPos = [1,3,2];
integrase = [1,0,0];
replication = [1,1,1];
death = [1,0,0];

data(:,1) = time;
data(:,2) = genFirstPos;
data(:,3) = integrase;
data(:,4) = replication;
data(:,5) = death;
S = 100; % size of a circle



for i = 1:length(time)
  
 
 %check if the cell did not die
    if death(i) == 0
        continue;
    else
        scatter(time(i),1,'x');
        break
    end
    
    % colour the circle based on which gene is in the first position
    hold on

 if genFirstPos(i) > 2. 
    scatter(time(i),1,S,'MarkerFaceColor','red');
 elseif genFirstPos(i) < 3 && genFirstPos(i)>1
    scatter(time(i),1,S,'MarkerFaceColor','blue');
 else
    scatter(time(i),1,S,'MarkerFaceColor','green');
 end
 
     if integrase (i) == 1
         plot([time(i),time(i+1)],[1,1],'b')
     else
         plot([time(i),time(i+1)],[1,1],'b--')
     end
 
end


