clc
clearvars
set(gcf,'units','points','position',[100,100,800,800])

%User Input Data
configuration = 2; %1=Linear, 2=Polygon
nozzles = 4;
side = 20;
rotation_speed = 2400;
print_speed = 1000;
ztravel_speed = 600;

%Initialise variables
print_time = 0;
xp = 0;
yp = 0;
zp = 0;
layer = -1;
z_val = 0;
m=1;
k=1;
u = 0;
angle= 0;


%Initialise print bed variables
rb = 100;
tb = 0:pi/180:1.8*pi;
xb = rb * cos(tb);
yb = rb * sin(tb);
zb = 0 * cos(tb);


if configuration==1
    ext_ang = pi;
else
    ext_ang = (2*pi)/nozzles;
    poly_rot = atan((sin(ext_ang/2)^2 - (nozzles-1)*(sin((nozzles-1)*ext_ang/2)^2))/(cos(ext_ang/2)*sin(ext_ang/2) - (nozzles-1)*cos((nozzles-1)*ext_ang/2)*sin((nozzles-1)*ext_ang/2)));
end


%Calculate Position of nozzles
x = nan(2*nozzles,2*nozzles);
y = nan(2*nozzles,2*nozzles);
z = nan(2*nozzles,2*nozzles);
tc = zeros(3, nozzles+1);
    for n=2:nozzles
        tc(1,n) = tc(1,n-1) + side*cos((n-1)*ext_ang);
        tc(2,n) = tc(2,n-1) + side*sin((n-1)*ext_ang);
    end

%Read the GCode file
input = fopen('Overlap3.gcode', 'r');
gcode = textscan(input,'%s','Delimiter','\n');
fclose(input);


hold on

%Define plot colors for each of the nozzles
Color = ["blue","green","red","cyan"];


for n = 1:nozzles
    p(m) = plot3(NaN,NaN,NaN,'Color',Color(n),'LineWidth',2);
    m=m+1;
end


%Plot the print bed, nozzles and the estimated print time annotation
pb= fill3(xb,yb,zb,[0.85 0.85 0.85]);
ht = plot3(tc(1,:), tc(2,:), tc(3,:),'LineWidth',5, 'Color', 'black','YDataSource','yt','XDataSource','xt','ZDataSource','zt');
axis([-110 110 -110 110])
txt = ['Estimated Print Time: ' num2str(round(print_time,2)) ' seconds'];
pt = text(30,100,txt,'FontSize',14);


%Parse the GCode
numCommands = numel(gcode{:});
for j = 1:numCommands
    if strncmp(gcode{1}{j},'G1 U',4)
        line = sscanf(gcode{1}{j}, ...
            '%c %f %c %f %c %f %c %f %c %f');

        values = zeros(1,5);
        values(1) = line(2);

        for i = 3:length(line)
            if strcmp(char(line(i)), 'U')
               values(2) = line(i+1);
            elseif strcmp(char(line(i)), 'X')
               values(3) = line(i+1);
            elseif strcmp(char(line(i)), 'Y')
               values(4) = line(i+1);
            elseif strcmp(char(line(i)), 'F')
               values(5) = line(i+1);
            end
        end

        %Update print_time with the current move and update variables
        print_dist = sqrt((values(3)-xp).^2 + (values(4)-yp).^2);
        print_time = print_time + (print_dist/(rotation_speed/60));
        xp = values(3);
        yp = values(4);

        %Rotate the plot based on 'U' axis value
        ang = values(2) - u;
        for r = 1:numel(p)
            rotate(p(r),[0 0 1],ang)
        end 
        rotate(pb,[0 0 1], ang)
        u = values(2);
       
       
        x = nan(2*nozzles,2*nozzles);
        y = nan(2*nozzles,2*nozzles);
        z = nan(2*nozzles,2*nozzles);
        axis([-110 110 -110 110])
        k=1;

        %Calculate the position of all the nozzles for the current move and update the plot
        a = zeros(1,nozzles);
        if configuration==1
           for n = 2:nozzles
               a(n) = a(n-1) + side*cos(ext_ang);
           end
        else
           for n = 2:nozzles
               a(n) = a(n-1) + side*cos((n-1)*ext_ang + pi/2 - angle);
           end
        end
        e_val = zeros(1,nozzles);
        for n = 1:nozzles
                if configuration==1
                    s = (values(3) + (by(1)-by(n))*cos(ext_ang) + (a(1)-a(n))*cos(ext_ang));
                    t = (values(4) + (by(1)-by(n))*sin(ext_ang) + (a(1)-a(n))*sin(ext_ang));
                else
                    s = (values(3) + by(n)*cos(angle) + (a(1)-a(n))*cos(pi/2 + angle));
                    t = (values(4) + by(n)*sin(angle) + (a(1)-a(n))*sin(pi/2 + angle));
                end
                x(n,k) = s;
                y(n,k) = t;
                z(n,k) = z_val;
                xt = tc(1,:) + values(3);    %t1+s;
                yt = tc(2,:) + values(4);     %t2+t;
                zt = tc(3,:);
                refreshdata(ht,'caller');
                drawnow;          
        end
        k=k+1;
      
    elseif strncmp(gcode{1}{j},'G1 X',4)
       
        line = sscanf(gcode{1}{j}, ...
            '%c %f %c %f %c %f %c %f %c %f %c %f %c %f');

        %define the array 'values'  
        values = zeros(1,7);
        values(1) = line(2);
        
        f=0;
        for i = 3:length(line)
            if strcmp(char(line(i)), 'X')
               values(2) = line(i+1);
            elseif strcmp(char(line(i)), 'Y')
               values(3) = line(i+1);
            elseif strcmp(char(line(i)), 'Z')
               values(4) = line(i+1);
            elseif strcmp(char(line(i)), 'E')
               values(5) = line(i+1);
            elseif strcmp(char(line(i)), ':')
               values(6+f) = line(i+1);
               f = f+1;
            end
        end

        %Update print_time with the current move and update variables
        print_dist = sqrt((values(2)-xp).^2 + (values(3)-yp).^2);
        print_time = print_time + (print_dist/(print_speed/60));
        xp = values(2);
        yp = values(3);
        
        %Update extrusion variables for each nozzle
        e_val = zeros(1,nozzles);
        for n = 1:nozzles
           e_val(n) = values(4+n);
        end
       
        %Calculate the position of all the nozzles for the current move and update the plot
        a = zeros(1,nozzles);
        if configuration==1
            for n = 2:nozzles
                a(n) = a(n-1) + side*cos(ext_ang);
            end
        else
            for n = 2:nozzles
                a(n) = a(n-1) + side*cos((n-1)*ext_ang + pi/2 - angle);
            end
        end
        for n = 1:nozzles
            if configuration==1
                    s = (values(2) + (by(1)-by(n))*cos(ext_ang) + (a(1)-a(n))*cos(ext_ang));
                    t = (values(3) + (by(1)-by(n))*sin(ext_ang) + (a(1)-a(n))*sin(ext_ang));
            else
                    s = (values(2) + by(n)*cos(angle) + (a(1)-a(n))*cos(pi/2 + angle));
                    t = (values(3) + by(n)*sin(angle) + (a(1)-a(n))*sin(pi/2 + angle));
            end
            x(n,k) = s;
            y(n,k) = t;
            z(n,k) = z_val;
            if e_val(n) > 0
                p(m) = plot3(x(n, k-1:k),y(n, k-1:k), z(n, k-1:k),'Color',Color(n),'LineWidth',2);
                m=m+1;
            end
            xt = tc(1,:) + values(2); %t1+s;
            yt = tc(2,:) + values(3); %t2+t;
            zt = tc(3,:);
            refreshdata(ht,'caller');
            drawnow;
        end
        k=k+1;
    elseif strncmp(gcode{1}{j},'G1 Z',4)
        %Increment the layer
        layer = layer + 1;

        line = sscanf(gcode{1}{j}, ...
            '%c %f %c %f %c %f %c %f %c %f %c %f %c %f');
        values = zeros(1,7);
        values(1) = line(2);
        z_val = line(4);
        tc(3,:) = tc(3,:) + line(4);
        z_dist = z_val - zp;

        %Update print_time with the current move and update variables
        print_time = print_time + (z_dist/(ztravel_speed/60));
        zp = z_val;
        

        %Update the position of the nozzles
        by = zeros(1,nozzles);
        if configuration==1
           for n = 2:nozzles
               by(n) = by(n-1) + side*sin(ext_ang);
           end
        else
           angle = poly_rot + layer*(pi-ext_ang);
           if angle<0
               angle = pi + angle;
           end
           for n = 2:nozzles
               by(n) = by(n-1) + side*sin((n-1)*ext_ang + pi/2 - angle);
           end
        end

    end

    %Update the print_time annotation
    delete(pt);
    txt = ['Estimated Print Time: ' num2str(round(print_time,2)) ' seconds'];
    pt = text(30,100,txt,'FontSize',14);
    pause(0.1)
end

hold off
