function f = filippova(f,used_N,wall_rotation,cbody,U,V)
%flippova 2nd order accurate treatments for a curved boundary.
% Idea: build unknown post collision distribution functions from near
% nodes proposed by Filippova et. al.
% Issue: unstable with high Reynolds number.

nx = Parameters.nx;
ny = Parameters.ny;
R  = Parameters.R;
Tau = Parameters.Tau;
w  = Parameters.w;
ex = Parameters.ex;
ey = Parameters.ey;

delta = zeros(nx,ny);
fstar = zeros(nx,ny,9);
UW    = zeros(nx,ny);
VW    = zeros(nx,ny);
Ubf   = zeros(nx,ny);
Vbf   = zeros(nx,ny);

for i=1:nx
   for j=1:ny
      if ~used_N(i,j)
         %% b=2
         if used_N(i,j+1)
            b=2;
            a=4;
            
            y= cbody(2) + sqrt( R^2 - (i-cbody(1))^2 );
            delta(i,j,a)= abs(j+1 - y);
            
            %%%%%%%%%%%collision%%%%%%%%%%%%%
            UW(i,j)= R * wall_rotation;
            VW(i,j)=0;
            
            if delta(i,j,a) >= 1/2
               Ubf(i,j)=(delta(i,j,a)-1)* U(i,j+1) / delta(i,j,a) + UW(i,j)/delta(i,j,a);
               Vbf(i,j)=(delta(i,j,a)-1)* V(i,j+1) / delta(i,j,a) + VW(i,j)/delta(i,j,a);
               X= (2*delta(i,j,a)-1)/Tau;
            end
            if delta(i,j,a) < 1/2
               Ubf(i,j)=U(i,j+1);
               Vbf(i,j)=V(i,j+1);
               X = (2*delta(i,j,a)-1)/(Tau - 1);
            end
            
            fstar(i,j,a) = w(a) * (1 + 3*(ex(a)*Ubf(i,j) + ey(a) * Vbf(i,j)) +   9/2*(   ex(a)*U(i,j+1)   +   ey(a) * V(i,j+1)  )^2 * (U(i,j+1)^2 + V(i,j+1)^2));
            
            f(i,j,b)= (1-X)*f(i,j+1,a) + X * fstar(i,j,a) + 2 * w(a)* 3 * (ex(b)*UW(i,j) + ey(b) * VW(i,j));
            %%%%%%%%%%%%%%%%%%%%%%%%%
         end
         %% 4
         if used_N(i,j-1)
            b=4;
            a=2;
            y= cbody(2) - sqrt( R^2 - (i-cbody(1))^2 );
            delta(i,j,a)= abs(y - (j-1));
            
            %%%%%%%%%%%collision%%%%%%%%%%%%%
            UW(i,j)= - R * wall_rotation;
            VW(i,j)=0;
            
            if delta(i,j,a) >= 1/2
               Ubf(i,j)=(delta(i,j,a)-1)* U(i,j-1) / delta(i,j,a) + UW(i,j)/delta(i,j,a);
               Vbf(i,j)=(delta(i,j,a)-1)* V(i,j-1) / delta(i,j,a) + VW(i,j)/delta(i,j,a);
               X= (2*delta(i,j,a)-1)/Tau;
            end
            if delta(i,j,a) < 1/2
               Ubf(i,j)=U(i,j-1);
               Vbf(i,j)=V(i,j-1);
               X = (2*delta(i,j,a)-1)/(Tau - 1);
            end
            
            fstar(i,j,a) = w(a) * (1 + 3*(ex(a)*Ubf(i,j) + ey(a) * Vbf(i,j)) + 9/2*(ex(a)*U(i,j-1) + ey(a) * V(i,j-1))^2 * (U(i,j-1)^2 + V(i,j-1)^2));
            
            f(i,j,b)= (1-X)*f(i,j-1,a) + X * fstar(i,j,a) + 2 * w(a)* 3 * (ex(b)*UW(i,j) + ey(b) * VW(i,j));
            %%%%%%%%%%%%%%%%%%%%%%%%%
         end
         %% 1
         if used_N(i+1,j)
            a=3;
            b=1;
            x= cbody(1) + sqrt( R^2 - (j-cbody(2))^2 ) ;
            delta(i,j,a)= abs(i+1 - x);
            
            %%%%%%%%%%%collision%%%%%%%%%%%%%
            UW(i,j)= 0;
            VW(i,j)= - R * wall_rotation;
            
            if delta(i,j,a) >= 1/2
               Ubf(i,j)=(delta(i,j,a)-1)* U(i+1,j) / delta(i,j,a) + UW(i,j)/delta(i,j,a);
               Vbf(i,j)=(delta(i,j,a)-1)* V(i+1,j) / delta(i,j,a) + VW(i,j)/delta(i,j,a);
               X= (2*delta(i,j,a)-1)/Tau;
            end
            if delta(i,j,a) < 1/2
               Ubf(i,j)=U(i+1,j);
               Vbf(i,j)=V(i+1,j);
               X = (2*delta(i,j,a)-1)/(Tau - 1);
            end
            
            fstar(i,j,a) = w(a) * (1 + 3*(ex(a)*Ubf(i,j) + ey(a) * Vbf(i,j)) + 9/2*(ex(a)*U(i+1,j) + ey(a) * V(i+1,j))^2 * (U(i+1,j)^2 + V(i+1,j)^2));
            
            f(i,j,b)= (1-X)*f(i+1,j,a) + X * fstar(i,j,a) + 2 * w(a)* 3 * (ex(b)*UW(i,j) + ey(b) * VW(i,j));
            %%%%%%%%%%%%%%%%%%%%%%%%%
         end
         %% 3
         if used_N(i-1,j)
            a=1;
            b=3;
            x= cbody(1) - sqrt( R^2 - (j-cbody(2))^2 ) ;
            delta(i,j,a)= abs(x - (i-1));
            
            %%%%%%%%%%%collision%%%%%%%%%%%%%
            UW(i,j)= 0;
            VW(i,j)= R * wall_rotation;
            
            if delta(i,j,a) >= 1/2
               Ubf(i,j)=(delta(i,j,a)-1)* U(i-1,j) / delta(i,j,a) + UW(i,j)/delta(i,j,a);
               Vbf(i,j)=(delta(i,j,a)-1)* V(i-1,j) / delta(i,j,a) + VW(i,j)/delta(i,j,a);
               X= (2*delta(i,j,a)-1)/Tau;
            end
            if delta(i,j,a) < 1/2
               Ubf(i,j)=U(i-1,j);
               Vbf(i,j)=V(i-1,j);
               X = (2*delta(i,j,a)-1)/(Tau - 1);
            end
            
            fstar(i,j,a) = w(a) * (1 + 3*(ex(a)*Ubf(i,j) + ey(a) * Vbf(i,j)) + 9/2*(ex(a)*U(i-1,j) + ey(a) * V(i-1,j))^2 * (U(i-1,j)^2 + V(i-1,j)^2));
            
            f(i,j,b)= (1-X)*f(i-1,j,a) + X * fstar(i,j,a) + 2 * w(a)* 3 * (ex(b)*UW(i,j) + ey(b) * VW(i,j));
            %%%%%%%%%%%%%%%%%%%%%%%%%
         end
         %% 5
         if used_N(i+1,j+1)
            b =   j-cbody(1)-cbody(2)-i;
            c =   1/2*( cbody(1)^2 + (j-cbody(2)-i)^2 - R^2 );
            x =   1/2*( -b + sqrt( b^2 - 4*c ) );
            y =   j + (x-i);
            a=7;
            bb=5;
            delta(i,j,a) = sqrt(( (i+1 - x)^2 + (j+1 - y)^2)  )/sqrt(2);
            
            %%%%%%%%%%%collision%%%%%%%%%%%%%
            teta= atan( abs(     (y - cbody(2)) / (x - cbody(1) )    )); % angel of the intersection of wall point
            UW(i,j)= abs( R * wall_rotation * cos(teta));
            VW(i,j)= -abs(R * wall_rotation * sin(teta));
            
            if delta(i,j,a) >= 1/2
               Ubf(i,j)=(delta(i,j,a)-1)* U(i+1,j+1) / delta(i,j,a) + UW(i,j)/delta(i,j,a);
               Vbf(i,j)=(delta(i,j,a)-1)* V(i+1,j+1) / delta(i,j,a) + VW(i,j)/delta(i,j,a);
               X= (2*delta(i,j,a)-1)/Tau;
            end
            if delta(i,j,a) < 1/2
               Ubf(i,j)=U(i+1,j+1);
               Vbf(i,j)=V(i+1,j+1);
               X = (2*delta(i,j,a)-1)/(Tau - 1);
            end
            
            fstar(i,j,a) = w(a) * (1 + 3*(ex(a)*Ubf(i,j) + ey(a) * Vbf(i,j)) + 9/2*(ex(a)*U(i+1,j+1) + ey(a) * V(i+1,j+1))^2 * (U(i+1,j+1)^2 + V(i+1,j+1)^2));
            
            f(i,j,bb)= (1-X)*f(i+1,j+1,a) + X * fstar(i,j,a) + 2 * w(a)* 3 * (ex(bb)*UW(i,j) + ey(bb) * VW(i,j));
            %%%%%%%%%%%%%%%%%%%%%%%%%
         end
         %% 7
         if used_N(i-1,j-1)
            b =   j-cbody(1)-cbody(2)-i;
            c =   1/2*( cbody(1)^2 + (j-cbody(2)-i)^2 - R^2 );
            x =   1/2*( -b - sqrt( b^2 - 4*c ) );
            y =   j + (x-i);
            a=5;
            bb=7;
            delta(i,j,a) = sqrt(1/2*( (x - (i-1) )^2 + (y - (j-1) )^2) );
            
            %%%%%%%%%%%collision%%%%%%%%%%%%%
            teta= atan( abs(     (y - cbody(2)) / (x - cbody(1) )    ));
            UW(i,j)= -abs(R * wall_rotation * cos(teta));
            VW(i,j)= abs(R * wall_rotation * sin(teta));
            
            if delta(i,j,a) >= 1/2
               Ubf(i,j)=(delta(i,j,a)-1)* U(i-1,j-1) / delta(i,j,a) + UW(i,j)/delta(i,j,a);
               Vbf(i,j)=(delta(i,j,a)-1)* V(i-1,j-1) / delta(i,j,a) + VW(i,j)/delta(i,j,a);
               X= (2*delta(i,j,a)-1)/Tau;
            end
            if delta(i,j,a) < 1/2
               Ubf(i,j)=U(i-1,j-1);
               Vbf(i,j)=V(i-1,j-1);
               X = (2*delta(i,j,a)-1)/(Tau - 1);
            end
            
            fstar(i,j,a) = w(a) * (1 + 3*(ex(a)*Ubf(i,j) + ey(a) * Vbf(i,j)) + 9/2*(ex(a)*U(i-1,j-1) + ey(a) * V(i-1,j-1))^2 * (U(i-1,j-1)^2 + V(i-1,j-1)^2));
            
            f(i,j,bb)= (1-X)*f(i-1,j-1,a) + X * fstar(i,j,a) + 2 * w(a)* 3 * (ex(bb)*UW(i,j) + ey(bb) * VW(i,j));
            %%%%%%%%%%%%%%%%%%%%%%%%%
         end
         %% 6
         if used_N(i-1,j+1)
            b =   -j-cbody(1)+cbody(2)-i;
            c =   1/2*(cbody(1)^2 +  (j-cbody(2)+i)^2 - R^2 );
            x =   1/2*( -b - sqrt( b^2 - 4*c ) );
            y =   j - (x-i);
            a=8;
            bb=6;
            delta(i,j,a) = sqrt(1/2*( (x - (i-1) )^2 + (j+1 -y)^2) );
            
            %%%%%%%%%%%collision%%%%%%%%%%%%%
            teta= atan( abs(     (y - cbody(2)) / (x - cbody(1) )    ));
            UW(i,j)= abs(R * wall_rotation * cos(teta));
            VW(i,j)= abs(R * wall_rotation * sin(teta));
            
            if delta(i,j,a) >= 1/2
               Ubf(i,j)=(delta(i,j,a)-1)* U(i-1,j+1) / delta(i,j,a) + UW(i,j)/delta(i,j,a);
               Vbf(i,j)=(delta(i,j,a)-1)* V(i-1,j+1) / delta(i,j,a) + VW(i,j)/delta(i,j,a);
               X= (2*delta(i,j,a)-1)/Tau;
            end
            if delta(i,j,a) < 1/2
               Ubf(i,j)=U(i-1,j+1);
               Vbf(i,j)=V(i-1,j+1);
               X = (2*delta(i,j,a)-1)/(Tau - 1);
            end
            
            fstar(i,j,a) = w(a) * (1 + 3*(ex(a)*Ubf(i,j) + ey(a) * Vbf(i,j)) + 9/2*(ex(a)*U(i-1,j+1) + ey(a) * V(i-1,j+1))^2 * (U(i-1,j+1)^2 + V(i-1,j+1)^2));
            
            f(i,j,bb)= (1-X)*f(i-1,j+1,a) + X * fstar(i,j,a) + 2 * w(a)* 3 * (ex(bb)*UW(i,j) + ey(bb) * VW(i,j));
            %%%%%%%%%%%%%%%%%%%%%%%%%
         end
         %% 8
         if used_N(i+1,j-1)
            b =   -j-cbody(1)+cbody(2)-i;
            c =   1/2*( cbody(1)^2 + (j-cbody(2)+i)^2 - R^2 );
            x =   1/2*( -b + sqrt( b^2 - 4*c ) );
            y =   j - (x-i);
            a=6;
            bb=8;
            delta(i,j,a) = sqrt(1/2*( (i+1 - x)^2 + (y - (j-1))^2) );
            
            %%%%%%%%%%%collision%%%%%%%%%%%%%
            teta= atan( abs(     (y - cbody(2)) / (x - cbody(1) )    ));
            UW(i,j)= -abs(R * wall_rotation * cos(teta));
            VW(i,j)= -abs(R * wall_rotation * sin(teta));
            
            if delta(i,j,a) >= 1/2
               Ubf(i,j)=(delta(i,j,a)-1)* U(i+1,j-1) / delta(i,j,a) + UW(i,j)/delta(i,j,a);
               Vbf(i,j)=(delta(i,j,a)-1)* V(i+1,j-1) / delta(i,j,a) + VW(i,j)/delta(i,j,a);
               X= (2*delta(i,j,a)-1)/Tau;
            end
            if delta(i,j,a) < 1/2
               Ubf(i,j)=U(i+1,j-1);
               Vbf(i,j)=V(i+1,j-1);
               X = (2*delta(i,j,a)-1)/(Tau - 1);
            end
            
            fstar(i,j,a) = w(a) * (1 + 3*(ex(a)*Ubf(i,j) + ey(a) * Vbf(i,j)) + 9/2*(ex(a)*U(i+1,j-1) + ey(a) * V(i+1,j-1))^2 * (U(i+1,j-1)^2 + V(i+1,j-1)^2));
            
            f(i,j,bb)= (1-X)*f(i+1,j-1,a) + X * fstar(i,j,a) + 2 * w(a)* 3 * (ex(bb)*UW(i,j) + ey(bb) * VW(i,j)) ;
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
         end
         
      end
   end
end

end