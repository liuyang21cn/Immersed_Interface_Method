    clear all

    icontour=1;
    iquiver=0;

    load xc.dat
    load yc.dat
    load surface1.dat
    load xsys.dat
    if iquiver==1 
      load uct.dat
      load vct.dat
    end
    if icontour==1
      load wot.dat
    end

    ms=1;    
    [n,m]=size(wot);
    nc=size(yc,1);
    mc=size(xc,1);
    ns=size(surface1,1);
    kn=n/nc;

% make_avi_movie_example2.m
% This is a simple example program to create an
% Audio Video Interleaved (AVI) movie that can be played
% independently of Matlab - see make_avi_movie_example1. 
% The movie is just a pseudocolor plot of a of a regular
% 2-D surface that is translated to make an animation. The
% example shows how to make a 2-D movie using the pcolor
% function in Matlab. Replacing pcolor with surf doesn't
% seem to work very well.

% Control no. of frames 
numframes = kn;

% Control how fast final movie will play; total time
% of movie in seconds will be (numframes)/(num_frames_per_second)
num_frames_per_second = 2;
dur = numframes/num_frames_per_second;
disp('                                              ');
disp('                                              ');
disp(['     duration of movie will be  ',num2str(dur), '  secs']);

% Create a new Audio Video Interleaved (AVI) file to be called
% something.avi using the avifile function
aviobj = avifile('flapper2d.avi','fps',num_frames_per_second); 

% Create plots of 2_D and add the 'frames' to the avi movie
% using the addframe function

for k = 1:numframes
      ks=(k-1)*nc+1;
      ke=ks+nc-1;
      H=pcolor(xc,yc,wot(ks:ke,:));
      %set(H,'EraseMode','xor');
      shading interp; 
      caxis([-4 4]);
      axis equal;
      hold on

      ks=(k-1)*ns+1;
      ke=ks+ns-1;
      for l=1:ms
        if l==1
          plot(xsys(ks:ke,1),xsys(ks:ke,2),'k-')
        else
          plot(xsys(ks:ke,2*l-1),xsys(ks:ke,2*l),'k-')
        end
        hold on
      end
      hold off

      frame = getframe(gca);
      aviobj = addframe(aviobj,frame);
end

% Tells matlab that there are no more frames to be added
% and to close the 2Dsinwave.avi file
aviobj = close(aviobj);
