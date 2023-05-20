function convergenceNeumann1
    hold on
    loglog(arrayfun(@(x) x^2,1:10));
    for w  = [3/2,3]
        drawThing = [];
        for  n = [4,9,19,39,79,99,119]
      %  for n = [4,9,19,39,49]
            [u_h,err] = solveNeumann(n,w,"false");
            drawThing = [drawThing,err]
        end
        loglog(drawThing);
    end
    loglog(arrayfun(@(x) x^2,drawThing(1):drawThing(1)+10));
end
