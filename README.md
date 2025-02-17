# Basic_Raycaster
Basic Raycaster Program

Build instructions: Compiled from the terminal, no makefile.
Compile: g++ basic_raycaster.cpp -o br
Run: ./br example_test_1.txt 

Example text file input: 

imsize 1024 512
eye 0 0 -10
viewdir 0 0 1
vfov 30
updir 0 1 0
bkgcolor 0.2 0.2 0.2
mtlcolor 1 0 0
sphere 0 0 15 2

mtlcolor 0 1 0
sphere 4 4 15 1

mtlcolor 0 0 1
sphere 10 4 15 1
sphere -10 -4 15 1
sphere 10 -4 15 1
sphere -10 4 15 1

