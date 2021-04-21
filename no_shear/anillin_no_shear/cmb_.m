for jj=0:1:100 % 3000
    jj
    
    run(['./outpt/out' num2str(jj)]);
    
    % grph ;
    plt ;

    axis([0 250 0 250]);
    axis equal
    axis off
    
    print( ['./mv/hex_' num2str(jj) '.png'],'-r200','-dpng' );

    clf
end