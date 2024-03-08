Welcome to my snowfall repo!
Here you'll find a handful of codes I wrote to simulate light snowfall, like a snow flurry.
I had developed this code for a class project around Christmas 2023. 

The code uses a finite difference solver to approximate the 2D Navier-Stokes equations to model wind patterns. 
The velocities of the wind, alongside some added downward velocity, determines the location of a handful of snowflakes. 
The example wind patterns I currently have coded are not physical, due to the limited time I had to complete this project for my class.
Regardless, the code gives a cute animation, as can be seen in the example video "snowfall.mp4".

The code is written to be used with openMP.

FILES:
Makefile -- generic makefile
getCPU.h -- supportive header file for CPU timing
parseCommand.h -- supportive header file to parse information from the command line
ins2d.C -- (main function) solves Navier Stokes equations, tracks snowflake particles in space and time, and outputs information to a Matlab file
snowMovie.m -- Matlab file to read snowflake particle information from ins2d.C and output a video

xmastree.jpeg -- background picture used in Matlab movie
snowfall.mp4 -- example movie


