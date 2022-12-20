%Initialize values
G = 1; %Gravitational Constant
N = 5000; %# of Timesteps  
tf = 2*pi;%Final time
t = linspace(0,tf,N+1); %Array of times
dt = t(2)-t(1); %Change between timesteps
n = 3; %Number of bodies

%Initial masses [Star a, Star b, Planet]
m = [1, 1, 0.0001]; 


%Position matrices 
x = zeros(n,N);
y = zeros(n,N);
z = zeros(n,N);

%Velocity matrices
u = zeros(n,N);
v = zeros(n,N);
w = zeros(n,N);

%Distances between stars
a = 1;
b = 1;

%Initial positions
x(:,1)= [a*0.5, -b*0.5, 0.4]; 
y(:,1)= [0,0,0];
z(:,1)= [0,0,0];

%Radius of planet a to star 
dxp = x(3)-x(1); 
dyp = y(3)-y(1);
dzp = z(3)-z(1);
r1 = sqrt(dxp^2+dyp^2+dzp^2);
sp = sqrt(G*(m(1)+m(3))/r1); 

%Initial velocities
w(:,1)=[0, 0, sp]; %For planet initial V = (U,V,W) + (u,v,w)*sqrt(GM/r)
v_calc_a = (m(2)/(m(1)+m(2)))*(sqrt((a/b)*(2*G*(m(1)+m(2)))/(a+b))); %Initial speed for Star a
v_calc_b = sqrt((b/a)*(2*G*(m(1)+m(2)))/(b+a))/2; %Initial speed for Star b
v(:,1) = [v_calc_a,-v_calc_b, v_calc_a]; 

%Loop 

for k=1:N-1
    for i= 1:n
        %Prevent overwriting of data
        u(i, k+1) = u(i, k);
        v(i, k+1) = v(i, k);
        w(i, k+1) = w(i, k);
        for j = 1:n
            if(i ~= j) %Each planet does not affect itself 
                %Updated distance between planets
                dx = x(j,k)-x(i,k);
                dy = y(j,k)-y(i,k);
                dz = z(j,k)-z(i,k);
                r = sqrt(dx^2+dy^2+dz^2);
                
                %Update velocities
                u(i, k+1) = u(i, k+1)+dt*G*m(j)*dx/r^3;
                v(i, k+1) = v(i, k+1)+dt*G*m(j)*dy/r^3;
                w(i, k+1) = w(i, k+1)+dt*G*m(j)*dz/r^3;
                
            end
        end
    end
    for i=1:n
        %Update positions
        x(i, k+1)=x(i, k)+dt*u(i, k+1);
        y(i, k+1)=y(i, k)+dt*v(i, k+1);
        z(i, k+1)=z(i, k)+dt*w(i, k+1);
        
    end
    
    if mod(k,5)==0
    %Plot
        cla
        scatter3(x(1, k+1), y(1,k+1),z(1,k+1), 'filled')
        hold on;
        scatter3(x(2, k+1), y(2,k+1),z(2,k+1), 'filled')
        scatter3(x(3, k+1), y(3,k+1),z(3,k+1), 'filled')
        plot3(x(1,1:k), y(1, 1:k), z(1, 1:k))
        plot3(x(2,1:k), y(2, 1:k),z(2, 1:k))
        plot3(x(3,1:k), y(3, 1:k),z(3, 1:k))
        axis equal
        axis([-2 2 -2 2 -2 2])
        
        drawnow
    end
end