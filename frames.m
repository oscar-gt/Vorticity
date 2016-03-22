%  Takes a (9 x 4096) matrix
function M = frames(w, tspan)

w = w';

for i = 1:length(tspan);
    wmat = reshape(w(:,i), 64,64);       %  Reshaping each time step vector
    %set(0,'DefaultFigureVisible','on'); %  Supressing figure output
    figure(i),
    
    pcolor(abs(wmat)),                   %  Making figure
    
    shading interp,
    colorbar,
    
    %drawnow,
    M(i) = getframe;             %  Storing to make a movie later
end

end