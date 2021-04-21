% run('./outpt/out0');

EdgesConsidered = [];
ListOfCellEdges = {};

CellEdgeBndryNode1 = [];
CellEdgeBndryNode2 = [];

for ii=1:1:length(Edge1)
    % ii

    CurrentEdge = 0;
    
    % find an edge not yet assigned to some cell edge
    for i=1:1:length(Edge1)
        if( sum( find(i==EdgesConsidered) )==0 ) % not yet considered
            if( (sum(Edge1(i)==Edge1)+sum(Edge1(i)==Edge2) == 3) | (sum(Edge2(i)==Edge1)+sum(Edge2(i)==Edge2) == 3) ) % each adjacent node adjacent to at most one edge
                EdgesConsidered = [EdgesConsidered i];
                continue;
            else
                CurrentEdge = i;
                break;
            end
        end
    end
    
    if( CurrentEdge==0 )
        break;
    end
    
    EdgeIdx = CurrentEdge;
    NodeIdx = Edge2(EdgeIdx);
    EdgeIndexesFrw = [EdgeIdx];
    
    for i=1:1:length(BndryX)
        
        if( sum(NodeIdx==Edge1)+sum(NodeIdx==Edge2) == 3 )
            CellEdgeBndryNode1 = [CellEdgeBndryNode1 NodeIdx];
            break
        end
        
        idxs1 = find( NodeIdx==Edge1 ) ;
        idxs2 = find( NodeIdx==Edge2 ) ;
        
        idxs = [idxs1 idxs2];
        
        for j=1:1:length(idxs)
            if( ~(idxs(j)==EdgeIdx) )
                NextEdgeIdx = idxs(j);
                
                if( Edge1(NextEdgeIdx)==NodeIdx )
                    NextNodeIdx = Edge2(NextEdgeIdx);
                else
                    NextNodeIdx = Edge1(NextEdgeIdx);
                end
                
                break
            end
        end
        
        EdgeIdx = NextEdgeIdx;
        NodeIdx = NextNodeIdx;
        
        EdgeIndexesFrw = [EdgeIndexesFrw EdgeIdx];
        
    end
    
    EdgeIdx = CurrentEdge;
    NodeIdx = Edge1(EdgeIdx);
    EdgeIndexesBckw = [];
    
    for i=1:1:length(BndryX)
        
        if( sum(NodeIdx==Edge1)+sum(NodeIdx==Edge2) == 3 )
            CellEdgeBndryNode2 = [CellEdgeBndryNode2 NodeIdx];
            break
        end
        
        idxs1 = find( NodeIdx==Edge1 ) ;
        idxs2 = find( NodeIdx==Edge2 ) ;
        
        idxs = [idxs1 idxs2];
        
        for j=1:1:length(idxs)
            if( ~(idxs(j)==EdgeIdx) )
                NextEdgeIdx = idxs(j);
                
                if( Edge1(NextEdgeIdx)==NodeIdx )
                    NextNodeIdx = Edge2(NextEdgeIdx);
                else
                    NextNodeIdx = Edge1(NextEdgeIdx);
                end
                
                break
            end
        end
        
        EdgeIdx = NextEdgeIdx;
        NodeIdx = NextNodeIdx;
        
        EdgeIndexesBckw = [EdgeIdx EdgeIndexesBckw];
    end
    
    EdgeIndexes = [EdgeIndexesBckw EdgeIndexesFrw];
    EdgesConsidered = [EdgesConsidered EdgeIndexesBckw EdgeIndexesFrw];
    
    ListOfCellEdges{1,size(ListOfCellEdges,2)+1} = EdgeIndexes;
end

% EdgeIndexes = ListOfCellEdges{1,140};
% 
% for i=1:1:length(EdgeIndexes)
%     j = EdgeIndexes(i);
%     plot( [BndryX(Edge1(j)) BndryX(Edge2(j))], [BndryY(Edge1(j)) BndryY(Edge2(j))], 'c' );
% end

ListOfCellEdgesFound = [];
ListsOfNodes = {};

for ii=1:1:size(ListOfCellEdges,2)
    
    % ii
    
    NotAllFound = 0;
    
    for i=1:1:size(ListOfCellEdges,2)
        if( ~sum(i==ListOfCellEdgesFound) )
            CurrentEdge0 = i;
            NotAllFound = 1;
            break;
        end
    end
    
    if(~NotAllFound)
        break;
    end
    
    CurrentNode0 = CellEdgeBndryNode1(CurrentEdge0);
    
    % CurrentEdge0 = 1;
    % CurrentNode0 = CellEdgeBndryNode1(CurrentEdge0);
    
    [ CurrentEdge1 CurrentEdge2 CurrentNode1 CurrentNode2 ] = getAdjacent(CurrentEdge0, CurrentNode0, CellEdgeBndryNode1, CellEdgeBndryNode2);
    
    [ CurrentEdge11 CurrentEdge12 CurrentNode11 CurrentNode12 ] = getAdjacent(CurrentEdge1, CurrentNode1, CellEdgeBndryNode1, CellEdgeBndryNode2);
    [ CurrentEdge21 CurrentEdge22 CurrentNode21 CurrentNode22 ] = getAdjacent(CurrentEdge2, CurrentNode2, CellEdgeBndryNode1, CellEdgeBndryNode2);
    
    [ CurrentEdge111 CurrentEdge112 CurrentNode111 CurrentNode112 ] = getAdjacent(CurrentEdge11, CurrentNode11, CellEdgeBndryNode1, CellEdgeBndryNode2);
    [ CurrentEdge121 CurrentEdge122 CurrentNode121 CurrentNode122 ] = getAdjacent(CurrentEdge12, CurrentNode12, CellEdgeBndryNode1, CellEdgeBndryNode2);
    [ CurrentEdge211 CurrentEdge212 CurrentNode211 CurrentNode212 ] = getAdjacent(CurrentEdge21, CurrentNode21, CellEdgeBndryNode1, CellEdgeBndryNode2);
    [ CurrentEdge221 CurrentEdge222 CurrentNode221 CurrentNode222 ] = getAdjacent(CurrentEdge22, CurrentNode22, CellEdgeBndryNode1, CellEdgeBndryNode2);
    
    CellEdge1 = [];
    CellEdge2 = [];
    CellEdge3 = [];
    CellEdge4 = [];
    
    if( CellEdgeBndryNode2(CurrentEdge0)==CurrentNode111 )
        CellEdge1 = ListOfCellEdges{1,CurrentEdge0};
        CellEdge2 = ListOfCellEdges{1,CurrentEdge1};
        CellEdge3 = ListOfCellEdges{1,CurrentEdge11};
        CellEdge4 = ListOfCellEdges{1,CurrentEdge111};
        
        ListOfCellEdgesFound = [ListOfCellEdgesFound CurrentEdge0 CurrentEdge1 CurrentEdge11 CurrentEdge111];
    end
    
    if( CellEdgeBndryNode2(CurrentEdge0)==CurrentNode112 )
        CellEdge1 = ListOfCellEdges{1,CurrentEdge0};
        CellEdge2 = ListOfCellEdges{1,CurrentEdge1};
        CellEdge3 = ListOfCellEdges{1,CurrentEdge11};
        CellEdge4 = ListOfCellEdges{1,CurrentEdge112};
        
        ListOfCellEdgesFound = [ListOfCellEdgesFound CurrentEdge0 CurrentEdge1 CurrentEdge11 CurrentEdge112];
    end
    
    if( CellEdgeBndryNode2(CurrentEdge0)==CurrentNode121 )
        CellEdge1 = ListOfCellEdges{1,CurrentEdge0};
        CellEdge2 = ListOfCellEdges{1,CurrentEdge1};
        CellEdge3 = ListOfCellEdges{1,CurrentEdge12};
        CellEdge4 = ListOfCellEdges{1,CurrentEdge121};
        
        ListOfCellEdgesFound = [ListOfCellEdgesFound CurrentEdge0 CurrentEdge1 CurrentEdge12 CurrentEdge121];
    end
    
    if( CellEdgeBndryNode2(CurrentEdge0)==CurrentNode122 )
        CellEdge1 = ListOfCellEdges{1,CurrentEdge0};
        CellEdge2 = ListOfCellEdges{1,CurrentEdge1};
        CellEdge3 = ListOfCellEdges{1,CurrentEdge12};
        CellEdge4 = ListOfCellEdges{1,CurrentEdge122};
        
        ListOfCellEdgesFound = [ListOfCellEdgesFound CurrentEdge0 CurrentEdge1 CurrentEdge12 CurrentEdge122];
    end
    
    if( CellEdgeBndryNode2(CurrentEdge0)==CurrentNode211 )
        CellEdge1 = ListOfCellEdges{1,CurrentEdge0};
        CellEdge2 = ListOfCellEdges{1,CurrentEdge2};
        CellEdge3 = ListOfCellEdges{1,CurrentEdge21};
        CellEdge4 = ListOfCellEdges{1,CurrentEdge211};
        
        ListOfCellEdgesFound = [ListOfCellEdgesFound CurrentEdge0 CurrentEdge2 CurrentEdge21 CurrentEdge211];
    end
    
    if( CellEdgeBndryNode2(CurrentEdge0)==CurrentNode212 )
        CellEdge1 = ListOfCellEdges{1,CurrentEdge0};
        CellEdge2 = ListOfCellEdges{1,CurrentEdge2};
        CellEdge3 = ListOfCellEdges{1,CurrentEdge21};
        CellEdge4 = ListOfCellEdges{1,CurrentEdge212};
        
        ListOfCellEdgesFound = [ListOfCellEdgesFound CurrentEdge0 CurrentEdge2 CurrentEdge21 CurrentEdge212];
    end
    
    if( CellEdgeBndryNode2(CurrentEdge0)==CurrentNode221 )
        CellEdge1 = ListOfCellEdges{1,CurrentEdge0};
        CellEdge2 = ListOfCellEdges{1,CurrentEdge2};
        CellEdge3 = ListOfCellEdges{1,CurrentEdge22};
        CellEdge4 = ListOfCellEdges{1,CurrentEdge221};
        
        ListOfCellEdgesFound = [ListOfCellEdgesFound CurrentEdge0 CurrentEdge2 CurrentEdge22 CurrentEdge221];
    end
    
    if( CellEdgeBndryNode2(CurrentEdge0)==CurrentNode222 )
        CellEdge1 = ListOfCellEdges{1,CurrentEdge0};
        CellEdge2 = ListOfCellEdges{1,CurrentEdge2};
        CellEdge3 = ListOfCellEdges{1,CurrentEdge22};
        CellEdge4 = ListOfCellEdges{1,CurrentEdge222};
        
        ListOfCellEdgesFound = [ListOfCellEdgesFound CurrentEdge0 CurrentEdge2 CurrentEdge22 CurrentEdge222];
    end
    
    if( Edge1(CellEdge1(1))==CurrentNode0 | Edge2(CellEdge1(1))==CurrentNode0 )
        CellEdge1 = flip(CellEdge1);
    end
    
    TerminalNodes1 = [ Edge1( CellEdge1(length(CellEdge1)) ) Edge2( CellEdge1(length(CellEdge1)) ) ];
    TerminalNodes2 = [ Edge1( CellEdge2(1) ) Edge2( CellEdge2(1) ) ];
    
    if( ~do_intersect(TerminalNodes1,TerminalNodes2) )
        CellEdge2 = flip(CellEdge2);
    end
    
    TerminalNodes1 = [ Edge1( CellEdge2(length(CellEdge2)) ) Edge2( CellEdge2(length(CellEdge2)) ) ];
    TerminalNodes2 = [ Edge1( CellEdge3(1) ) Edge2( CellEdge3(1) ) ];
    
    if( ~do_intersect(TerminalNodes1,TerminalNodes2) )
        CellEdge3 = flip(CellEdge3);
    end
    
    TerminalNodes1 = [ Edge1( CellEdge3(length(CellEdge3)) ) Edge2( CellEdge3(length(CellEdge3)) ) ];
    TerminalNodes2 = [ Edge1( CellEdge4(1) ) Edge2( CellEdge4(1) ) ];
    
    if( ~do_intersect(TerminalNodes1,TerminalNodes2) )
        CellEdge4 = flip(CellEdge4);
    end
    
    ListOfNodes = [getNodes(Edge1(CellEdge1),Edge2(CellEdge1)) getNodes(Edge1(CellEdge2),Edge2(CellEdge2)) getNodes(Edge1(CellEdge3),Edge2(CellEdge3)) getNodes(Edge1(CellEdge4),Edge2(CellEdge4))] ;
    
    ListsOfNodes{1,size(ListsOfNodes,2)+1} = ListOfNodes;
    
    % patch( BndryX(ListOfNodes), BndryY(ListOfNodes), 'c','EdgeColor','none' );
    
end



BasalEdgeIdxs = find(EdgeId==3);

BasalNodes1 = Edge1( BasalEdgeIdxs );
BasalNodes2 = Edge2( BasalEdgeIdxs );

CurrentEdge = BasalEdgeIdxs(1);
CurrentNode = Edge1(CurrentEdge);

BasalEdgesReordered = [CurrentEdge];

for i=1:1:length(BasalNodes1)-1
    idxs = [find(BasalNodes1==CurrentNode) find(BasalNodes2==CurrentNode)];
    
    if( BasalEdgeIdxs(idxs(1))==CurrentEdge )
        NextEdge = BasalEdgeIdxs(idxs(2));
    else
        NextEdge = BasalEdgeIdxs(idxs(1));
    end
    
    if( Edge1(NextEdge)==CurrentNode )
        NextNode = Edge2(NextEdge);
    else
        NextNode = Edge1(NextEdge);
    end
    
    BasalEdgesReordered = [BasalEdgesReordered NextEdge];
    
    CurrentEdge = NextEdge;
    CurrentNode = NextNode;
end

BasalNodes = getNodes(Edge1(BasalEdgesReordered),Edge2(BasalEdgesReordered));

% patch(BndryX(BasalNodes),BndryY(BasalNodes),'b','EdgeColor','none')

hold on

for i=1:1:size(ListsOfNodes,2)
    patch( BndryX(ListsOfNodes{1,i}), BndryY(ListsOfNodes{1,i}), 'c','EdgeColor','none' );
end

% patch(BndryX(BasalNodes),BndryY(BasalNodes),'b','EdgeColor','none')
% patch(BndryX(BasalNodes),BndryY(BasalNodes),'c','EdgeColor','none')
patch('Faces',1:1:length(BndryX(BasalNodes)),'Vertices',[BndryX(BasalNodes)' BndryY(BasalNodes)'],'FaceColor',[0.65 1 1],'EdgeColor','none') ;