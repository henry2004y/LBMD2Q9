% Collection of unused codes

%%
   % Calculate macroscopic velocities
   for j=1:ny; for i=1:nx
      if used_N(i,j)
         Rhotemp = 0;
         MomtempX = 0;
         MomtempY = 0;
         for iDir=1:9
            Rhotemp = Rhotemp + f(i,j,iDir);
            MomtempX = MomtempX + f(i,j,iDir)*ex(iDir);
            MomtempY = MomtempY + f(i,j,iDir)*ey(iDir);
         end
        
         Rho(i,j) = Rhotemp;
         U(i,j) = MomtempX / Rho(i,j);
         V(i,j) = MomtempY / Rho(i,j);
      end  
   end;end

%%
   for j=1:ny; for i=1:nx
      if ~used_N(i,j)
         U(i,j) = 0;
         V(i,j) = 0;
      end
   end; end