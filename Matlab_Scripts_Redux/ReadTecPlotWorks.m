function data = ReadTecPlot(filename, PlaneName)
id=fopen(filename);
%PlaneName='Ground'
Space=20-length(PlaneName)
PlaneName=[PlaneName, blanks(Space)]
r=1;
t=1;
j=0;
%while 1
readin=(fgetl(id))
while 1
    if strcmp(readin,[' ZONE T="',PlaneName,'"'])
        for i=1:4
            if i==1
                j=j+1
            end
            readin1=fgetl(id);
        end
        i=1;
        while 1
            readin=fgetl(id);   
            if strncmp(readin,' ZONE',5)
                break
            end
            if readin==-1 
                break
            end
            data(i,:,t) = str2num(readin);
            i=i+1;
        end
        t=t+1;
    else
        readin=fgetl(id);
    end
    if readin==-1 
        break
    end
    r=r+1;
end


