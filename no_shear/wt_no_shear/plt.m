% hold on
% for i=1:1:size(BndryX,1)
%     plot(BndryX(i,:),BndryY(i,:),'.-k','LineWidth',2,'MarkerSize',12);
% end

% Edge1 = load('./meshfiles/edge_1.dat');
% Edge2 = load('./meshfiles/edge_2.dat');

% plot(BndryX,BndryY,'.k');

hold on
for i=1:1:length(Edge1)
    if( EdgeId(i)==3 )
        % continue;
    end
    
    % plot( [BndryX(Edge1(i)) BndryX(Edge2(i))], [BndryY(Edge1(i)) BndryY(Edge2(i))], 'r' );
    
    plot( [BndryX(Edge1(i)) BndryX(Edge2(i))], [BndryY(Edge1(i)) BndryY(Edge2(i))], '.-k' );
    
end
