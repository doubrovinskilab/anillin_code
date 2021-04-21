function [ CurrentEdge1 CurrentEdge2 CurrentNode1 CurrentNode2 ] = getAdjacent(CurrentEdge0, CurrentNode0, CellEdgeBndryNode1, CellEdgeBndryNode2)

% All (3) edges adjacent to node CorrentNode
CurrentNodeActual = CurrentNode0;
CurrentEdgeActual = CurrentEdge0;

AdjacentCellEdges = [find( CurrentNodeActual==CellEdgeBndryNode1 ) find( CurrentNodeActual==CellEdgeBndryNode2 )];

% (2) adjacent edges except CurrentEdge
AdjacentCellEdges = AdjacentCellEdges( find((AdjacentCellEdges==CurrentEdgeActual)-1) );

CurrentEdgeNext1 = AdjacentCellEdges(1);
CurrentEdgeNext2 = AdjacentCellEdges(2);

if( CellEdgeBndryNode1(CurrentEdgeNext1)==CurrentNodeActual )
    CurrentNodeNext1 = CellEdgeBndryNode2(CurrentEdgeNext1);
else
    CurrentNodeNext1 = CellEdgeBndryNode1(CurrentEdgeNext1);
end

if( CellEdgeBndryNode1(CurrentEdgeNext2)==CurrentNodeActual )
    CurrentNodeNext2 = CellEdgeBndryNode2(CurrentEdgeNext2);
else
    CurrentNodeNext2 = CellEdgeBndryNode1(CurrentEdgeNext2);
end

CurrentEdge1 = CurrentEdgeNext1 ;
CurrentEdge2 = CurrentEdgeNext2 ;

CurrentNode1 = CurrentNodeNext1 ;
CurrentNode2 = CurrentNodeNext2 ;

end