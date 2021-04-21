%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab code to generate the initial geometry of embryonic tissue 
% x.dat and y.dat contain the coordinates of the nodes discretizing 
% the profiles of cellular membranes
% edge_1.dat edge_2.dat describe edges: nodes edge_1(i) and edge_2(i)
% are adjacent to i:th edge
% apical_pulled.dat mesodermal apical edges (to be subjected to elevated
% stress)
% lateral_pulled.dat mesodermal lateral edges
% basal_edges.dat basal edges (all of them, mesodermal and not)
% apical_points.dat nodes on the basal boundary
% all output written to folder ./meshfiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NumPoints = 400*2 ;
NumCells = 100 ;

PointsPerCell = floor(NumPoints/NumCells);

R = 1.5*75.0; % 100.0;
L = 1.5*30.0 ;

xs = [];
ys = [];

for i=1:1:NumPoints
    Theta = 2*pi*(i-1)/NumPoints;
    
    xs = [xs R*cos(Theta)];
    ys = [ys R*sin(Theta)];
end

ApicalPoints = 1:1:NumPoints;

dx = xs(2)-xs(1);
dy = ys(2)-ys(1);
dl = sqrt(dx*dx+dy*dy);

NumPointsInRay = floor(L/dl);

for i=1:1:NumPoints
    Theta = 2*pi*(i-1)*PointsPerCell/NumPoints;

    if( Theta>=2*pi )
        break;
    end
    
    for j=1:1:NumPointsInRay
        xs = [xs (R - j*dl)*cos(Theta)];
        ys = [ys (R - j*dl)*sin(Theta)];
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
NumBottomPoints = 3*2;
NumCells = ( length(xs)-NumPoints )/NumPointsInRay ;
R1 = R - NumPointsInRay*dl ;
BottomPointsX = [];
BottomPointsY = [];

for i=1:1:NumBottomPoints*NumCells
    
    if ( (i-1)/NumBottomPoints-floor((i-1)/NumBottomPoints) == 0 )
        continue ;
    end
    
    Theta = 2*pi*(i-1)/(NumBottomPoints*NumCells);
    
    BottomPointsX = [BottomPointsX R1*cos(Theta)];
    BottomPointsY = [BottomPointsY R1*sin(Theta)];
end

dx = BottomPointsX(2)-BottomPointsX(1);
dy = BottomPointsY(2)-BottomPointsY(1);
dl_ = sqrt(dx*dx+dy*dy);

xs = [xs BottomPointsX];
ys = [ys BottomPointsY];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Edge1 = [];
Edge2 = [];

for i=1:1:length(xs)
    % for j=1:1:length(xs)
    for j=i:1:length(xs)
        
        rx = xs(j)-xs(i);
        ry = ys(j)-ys(i);
        r = sqrt( rx*rx + ry*ry );
        
        if( r<1.01*dl & r>0.99*dl )
            Edge1 = [Edge1 i];
            Edge2 = [Edge2 j];
        end
        
        if( r<1.01*dl_ & r>0.99*dl_ )
            Edge1 = [Edge1 i];
            Edge2 = [Edge2 j];
        end
    end
end

ApicalPulled = [];
LateralPulled = [];

for i=1:1:NumPoints
    % Theta = 2*pi*(i-1)/NumPoints;
    
    if( i>=PointsPerCell*17+1 )
        if( i<=PointsPerCell*33+1 )
            ApicalPulled = [ApicalPulled i];
        end
    end
end

for i=1:1:length(xs)
    
    if( i>=(NumPoints+NumPointsInRay*17+1) )
        if( i<=(NumPoints+NumPointsInRay*(33+1)) )
            LateralPulled = [LateralPulled i];
        end
    end
    
end

ApicalEdgesPulled = [];
LateralEdgesPulled = [];

for i=1:1:length(ApicalPulled)
    for j=1:1:length(ApicalPulled)
        rx = xs(ApicalPulled(j))-xs(ApicalPulled(i));
        ry = ys(ApicalPulled(j))-ys(ApicalPulled(i));
        r = sqrt( rx*rx+ry*ry );
        
        if( r<1.1*dl & r>0.9*dl )
            for k=1:1:length(Edge1)
                if( (Edge1(k)==ApicalPulled(i) & Edge2(k)==ApicalPulled(j)) | (Edge2(k)==ApicalPulled(i) & Edge1(k)==ApicalPulled(j)) )
                    ApicalEdgesPulled = [ApicalEdgesPulled k];
                end
            end
        end
    end
end

for i=1:1:length(LateralPulled)
    for j=1:1:length(LateralPulled)
        rx = xs(LateralPulled(j))-xs(LateralPulled(i));
        ry = ys(LateralPulled(j))-ys(LateralPulled(i));
        r = sqrt( rx*rx+ry*ry );
        
        if( r<1.1*dl & r>0.9*dl )
            for k=1:1:length(Edge1)
                if( (Edge1(k)==LateralPulled(i) & Edge2(k)==LateralPulled(j)) | (Edge2(k)==LateralPulled(i) & Edge1(k)==LateralPulled(j)) )
                    LateralEdgesPulled = [LateralEdgesPulled k];
                end
            end
        end
    end
end

for i=1:1:length(LateralPulled)
    for j=1:1:length(ApicalPulled)
        rx = xs(LateralPulled(i))-xs(ApicalPulled(j));
        ry = ys(LateralPulled(i))-ys(ApicalPulled(j));
        r = sqrt( rx*rx+ry*ry );
        
        if( r<1.1*dl & r>0.9*dl )
            for k=1:1:length(Edge1)
                if( (Edge1(k)==LateralPulled(i) & Edge2(k)==ApicalPulled(j)) | (Edge2(k)==LateralPulled(i) & Edge1(k)==ApicalPulled(j)) )
                    LateralEdgesPulled = [LateralEdgesPulled k];
                end
            end
        end
    end
end

BasalEdges = [];

for i=1:1:length(Edge1)
    rx1 = xs(Edge1(i));
    ry1 = ys(Edge1(i));
    
    rx2 = xs(Edge2(i));
    ry2 = ys(Edge2(i));
    
    r1 = sqrt(rx1*rx1 + ry1*ry1);
    r2 = sqrt(rx2*rx2 + ry2*ry2);
    
    if( r1<R1+0.1*dl & r2<R1+0.1*dl )
        BasalEdges = [BasalEdges i];
    end
end

xs = xs + 128; % R;
ys = ys + 128; % R;

plot(xs,ys,'.k');
hold on

for i=1:1:length(Edge1)
    idx1 = Edge1(i);
    idx2 = Edge2(i);
    
    plot([xs(idx2) xs(idx1)],[ys(idx2) ys(idx1)],'r');
end

for i=1:1:length(ApicalEdgesPulled)
    idx1 = Edge1(ApicalEdgesPulled(i));
    idx2 = Edge2(ApicalEdgesPulled(i));
    
    plot([xs(idx2) xs(idx1)],[ys(idx2) ys(idx1)],'c');
end

for i=1:1:length(LateralEdgesPulled)
    idx1 = Edge1(LateralEdgesPulled(i));
    idx2 = Edge2(LateralEdgesPulled(i));
    
    plot([xs(idx2) xs(idx1)],[ys(idx2) ys(idx1)],'b');
end

for i=1:1:length(BasalEdges)
    idx1 = Edge1(BasalEdges(i));
    idx2 = Edge2(BasalEdges(i));
    
    plot([xs(idx2) xs(idx1)],[ys(idx2) ys(idx1)],'c');
end

plot(xs(ApicalPulled),ys(ApicalPulled),'.g');
plot(xs(LateralPulled),ys(LateralPulled),'.m');

% Uncomment to write output:

x = xs';
y = ys';

save('./meshfiles/x.dat','x','-ascii');
save('./meshfiles/y.dat','y','-ascii');

Edge1 = Edge1';
Edge2 = Edge2';

save('./meshfiles/edge_1.dat','Edge1','-ascii');
save('./meshfiles/edge_2.dat','Edge2','-ascii');

ApicalEdgesPulled = ApicalEdgesPulled';
LateralEdgesPulled = LateralEdgesPulled';
BasalEdges = BasalEdges';
ApicalPoints = ApicalPoints';

save('./meshfiles/apical_pulled.dat','ApicalEdgesPulled','-ascii');
save('./meshfiles/lateral_pulled.dat','LateralEdgesPulled','-ascii');
save('./meshfiles/basal_edges.dat','BasalEdges','-ascii');
save('./meshfiles/apical_points.dat','ApicalPoints','-ascii');