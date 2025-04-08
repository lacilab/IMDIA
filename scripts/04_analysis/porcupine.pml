# Load the GROMACS .gro file (assuming it's in the current directory)
load extreme1.gro, extreme1

# Your original commands
hide ribbon
show cartoon
split_states extreme1, 1, 2
split_states extreme1, 2000, 2001
set_name extreme1_0001, 0001
set_name extreme1_2001, 2001
delete extreme1_0002
delete extreme1_2000
run modevectors.py
modevectors 0001, 2001, cutoff=0, head_length=1, head=0.4, headrgb=(1, 0.2, 0.1), tailrgb=(1, 0.2, 0.1), notail=0

# Visualization settings
set cartoon_ring_mode, 3
set ray_shadows,0
set ray_trace_mode,1
bg_color white
set antialias,2
set ray_trace_color, black
set ray_trace_gain, 0.005

set_view (\
     0.153344378,   -0.861105382,   -0.484749764,\
     0.879321575,   -0.104914077,    0.464531869,\
    -0.450869858,   -0.497483373,    0.741099834,\
     0.000018567,   -0.000001550, -886.166992188,\
    -6.509384632,    7.209011555,    4.623818398,\
   509.141479492, 1263.196044922,   -8.000000000 )

# Render and save image
ray 1920, 1080
png porcupine.png 

# Close PyMOL
quit
