

figure(1);
contourf( MCDS.mesh.X , MCDS.mesh.Y , MCDS.continuum_variables(1).data ) ; axis image; colorbar


figure(2); 
plot( MCDS.discrete_cells.custom.unbound_external_ACE2 + MCDS.discrete_cells.custom.bound_external_ACE2 + ...
    MCDS.discrete_cells.custom.unbound_internal_ACE2 + MCDS.discrete_cells.custom.bound_internal_ACE2  );
title( 'total receptor' ); 


figure(3);
plot( MCDS.discrete_cells.custom.virion ); 
title( 'virion' )

figure(4); 
plot( MCDS.discrete_cells.custom.uncoated_virion );
title( 'uncoated' )

figure(5); 
plot( MCDS.discrete_cells.custom.viral_RNA );
title( 'RNA'  )

figure(6); 
plot( MCDS.discrete_cells.custom.viral_protein );
title( 'protein'  )

figure(7); 
plot( MCDS.discrete_cells.custom.assembled_virion );
title( 'assembled virion'  )

figure(8); 
plot( MCDS.discrete_cells.custom.virion + MCDS.discrete_cells.custom.uncoated_virion + MCDS.discrete_cells.custom.viral_RNA );
title( 'virion + uncoated + RNA'  )

disp( 'min virion concentration' );
min(min( MCDS.continuum_variables(1).data ))