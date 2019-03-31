function f = bozidi(f,used_N,cbody)

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



for j=1:ny
   for i=1:nx
      if ~used_N(i,j)
         %% b=2
         if used_N(i,j+1) == 1
            b=2;
            a=4;
            
            y= cbody(2) + sqrt( R*R - (i-cbody(1))*(i-cbody(1)) ) ;
            delta(i,j,a)= abs(j+1 - y);
            
            %%%%%%%%%%%collision%%%%%%%%%%%%%
            if delta(i,j,a) >= 1/2
               f(i,j,b)= (1/(2*delta(i,j,a)))*f(i,j+1,a) + (2*delta(i,j,a)-1)/(2*delta(i,j,a))* f(i,j+1,b);
               
            end
            if delta(i,j,a) < 1/2
               f(i,j,b)= 2*delta(i,j,a)*  f(i,j+1,a) + ( 1 - 2*delta(i,j,a) )* f(i,j+2,b);
               
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
         end
         %% 4
         if used_N(i,j-1)
            b=4;
            a=2;
            y= cbody(2) - sqrt( R*R - (i-cbody(1))*(i-cbody(1)) ) ;
            delta(i,j,a)= abs(y - (j-1));
            
            %%%%%%%%%%%collision%%%%%%%%%%%%%
            if delta(i,j,a) >= 1/2
               f(i,j,b)= (1/(2*delta(i,j,a)))*f(i,j-1,a) + (2*delta(i,j,a)-1)/(2*delta(i,j,a))* f(i,j-1,b);
               
            end
            if delta(i,j,a) < 1/2
               f(i,j,b)= 2*delta(i,j,a)*  f(i,j-1,a) + ( 1 - 2*delta(i,j,a) )* f(i,j-2,b);
               
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%
         end
         %% 1
         if used_N(i+1,j)
            a=3;
            b=1;
            x= cbody(1) + sqrt( R*R - (j-cbody(2))*(j-cbody(2)) ) ;
            delta(i,j,a)= abs(i+1 - x);
            
            %%%%%%%%%%%collision%%%%%%%%%%%%%
            if delta(i,j,a) >= 1/2
               f(i,j,b)= (1/(2*delta(i,j,a)))*f(i+1,j,a) + (2*delta(i,j,a)-1)/(2*delta(i,j,a))* f(i+1,j,b);
               
            end
            if delta(i,j,a) < 1/2
               f(i,j,b)= 2*delta(i,j,a)*  f(i+1,j,a) + ( 1 - 2*delta(i,j,a) )* f(i+2,j,b);
               
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%
         end
         %% 3
         if used_N(i-1,j)
            a=1;
            b=3;
            x= cbody(1) - sqrt( R*R - (j-cbody(2))*(j-cbody(2)) ) ;
            delta(i,j,a)= abs(x - (i-1));
            
            %%%%%%%%%%%collision%%%%%%%%%%%%%
            if delta(i,j,a) >= 1/2
               f(i,j,b)= (1/(2*delta(i,j,a)))*f(i-1,j,a) + (2*delta(i,j,a)-1)/(2*delta(i,j,a))* f(i-1,j,b);
               
            end
            if delta(i,j,a) < 1/2
               f(i,j,b)= 2*delta(i,j,a)*  f(i-2,j,a) + ( 1 - 2*delta(i,j,a) )* f(i-2,j,b);
               
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%
         end
         %% 5
         if used_N(i+1,j+1)
            b =   j-cbody(1)-cbody(2)-i;
            c =   1/2*( cbody(1)*cbody(1) + (j-cbody(2)-i)*(j-cbody(2)-i) - R*R );
            x =   1/2*( -b + sqrt( b*b - 4*c ) );
            y =   j + (x-i);
            a=7;
            bb=5;
            delta(i,j,a) = sqrt(( (i+1 - x)*(i+1 - x) + (j+1 - y)*(j+1 - y))  )/sqrt(2);
            
            %%%%%%%%%%%collision%%%%%%%%%%%%%
            if delta(i,j,a) >= 1/2
               f(i,j,bb)= (1/(2*delta(i,j,a)))*f(i+1,j+1,a) + (2*delta(i,j,a)-1)/(2*delta(i,j,a))* f(i+1,j+1,bb);
               
            end
            if delta(i,j,a) < 1/2
               f(i,j,bb)= 2*delta(i,j,a)*  f(i+1,j+1,a) + ( 1 - 2*delta(i,j,a) )* f(i+2,j+2,bb);
               
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%
         end
         %% 7
         if used_N(i-1,j-1)
            b =   j-cbody(1)-cbody(2)-i;
            c =   1/2*( cbody(1)*cbody(1) + (j-cbody(2)-i)*(j-cbody(2)-i) - R*R );
            x =   1/2*( -b - sqrt( b*b - 4*c ) );
            y =   j + (x-i);
            a=5;
            bb=7;
            delta(i,j,a) = sqrt(1/2*( (x - (i-1) )*(x - (i-1) ) + (y - (j-1) )*(y - (j-1) )) );
            
            %%%%%%%%%%%collision%%%%%%%%%%%%%
            if delta(i,j,a) >= 1/2
               f(i,j,bb)= (1/(2*delta(i,j,a)))*f(i-1,j-1,a) + (2*delta(i,j,a)-1)/(2*delta(i,j,a))* f(i-1,j-1,bb);
               
            end
            if delta(i,j,a) < 1/2
               f(i,j,bb)= 2*delta(i,j,a)*  f(i-1,j-1,a) + ( 1 - 2*delta(i,j,a) )* f(i-2,j-2,bb);
               
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%
         end
         %% 6
         if used_N(i-1,j+1)
            b =   -j-cbody(1)+cbody(2)-i;
            c =   1/2*(cbody(1)*cbody(1) +  (j-cbody(2)+i)*(j-cbody(2)+i) - R*R );
            x =   1/2*( -b - sqrt( b*b - 4*c ) );
            y =   j - (x-i);
            a=8;
            bb=6;
            delta(i,j,a) = sqrt(1/2*( (x - (i-1) )*(x - (i-1) ) + (j+1 -y)*(j+1 -y)) );
            
            %%%%%%%%%%%collision%%%%%%%%%%%%%
            if delta(i,j,a) >= 1/2
               f(i,j,bb)= (1/(2*delta(i,j,a)))*f(i-1,j+1,a) + (2*delta(i,j,a)-1)/(2*delta(i,j,a))* f(i-1,j+1,bb);
               
            end
            if delta(i,j,a) < 1/2
               f(i,j,bb)= 2*delta(i,j,a)*  f(i-1,j+1,a) + ( 1 - 2*delta(i,j,a) )* f(i-2,j+2,bb);
               
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%
         end
         %% 8
         if used_N(i+1,j-1)
            b =   -j-cbody(1)+cbody(2)-i;
            c =   1/2*( cbody(1)*cbody(1) + (j-cbody(2)+i)*(j-cbody(2)+i) - R*R );
            x =   1/2*( -b + sqrt( b*b - 4*c ) );
            y =   j - (x-i);
            a=6;
            bb=8;
            delta(i,j,a) = sqrt(1/2*( (i+1 - x)*(i+1 - x) + (y - (j-1))*(y - (j-1))) );
            
            %%%%%%%%%%%collision%%%%%%%%%%%%%
            if delta(i,j,a) >= 1/2
               f(i,j,bb)= (1/(2*delta(i,j,a)))*f(i+1,j-1,a) + (2*delta(i,j,a)-1)/(2*delta(i,j,a))* f(i+1,j-1,bb);
               
            end
            if delta(i,j,a) < 1/2
               f(i,j,bb)= 2*delta(i,j,a)*  f(i+1,j-1,a) + ( 1 - 2*delta(i,j,a) )* f(i+2,j-2,bb);
               
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
         end
         
      end
   end
end
