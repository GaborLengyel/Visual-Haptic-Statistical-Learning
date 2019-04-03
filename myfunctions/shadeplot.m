function [Hp,Hs]=shadeplot(x,y,sd,linestyle,shadecolor)
% function [Hp,Hs]=shadeplot(x,y,sd,linestyle,shadecolor)
%
% Plot x,y with shaded regions to indicate y +/- sd.
%
% x,y,sd - Vectors of identical length.
% linestyle - Attributes as per standard plot function.
% shadecolor - RGB vector for color of shaded region.
%
% Returns handles to the plot (Hp) and shaded region (Hs).

%added by dmw to deal with nans
w=isnan(y)|isnan(sd);
x(w)=[];
y(w)=[];
sd(w)=[];

    k = 0;
    for i=1:1:length(x)
        k = k + 1;
        X(k) = x(i);
        Y(k) = y(i)-sd(i);
    end

    for i=length(x):-1:1
        k = k + 1;
        X(k) = x(i);
        Y(k) = y(i)+sd(i);
    end

    k = k + 1;
    X(k) = x(1);
    Y(k) = y(1)-sd(1);
    
    if( size(shadecolor,2) == 1 )
        shadecolor = shadecolor';
    end

    hold on;
    Hs = fill(X,Y,shadecolor,'LineStyle','none');
    Hp = plot(x,y,linestyle);

end
