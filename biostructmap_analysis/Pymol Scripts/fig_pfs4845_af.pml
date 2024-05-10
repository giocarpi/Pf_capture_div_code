## show AF confidence

# background white and alphachecker off
bg_color white
set ray_opaque_background, off

load ../biostructmap_input/pfs4845_AF_cleaned.pdb
load ../tajima_d_models/help_scale_AF.pdb

# set view (orient + zoom center, 60)
set_view (\
     0.904379010,    0.057848230,    0.422790945,\
    -0.420346409,   -0.049968023,    0.905986845,\
     0.073535763,   -0.997074127,   -0.020873735,\
     0.000000000,    0.000000000, -340.276916504,\
    -1.662139893,   -0.374973297,   -3.868373871,\
   268.276916504,  412.276916504,  -20.000000000 )

# run quick_color2
as cartoon

# turn off help_scale2
hide everything, help_scale_AF

#Select all residues with a mapped data value. Can change the default 'no-value'
#option when writing to pdb b factor using biostructmap if needed.

spectrum b

# Taking picture
unset specular
set ray_trace_gain, 1
set ray_trace_mode, 0
bg_color white
set ray_trace_color, black
unset depth_cue
set antialias, 2
set ambient, 0.2
set direct, 0.9

png ../Pymol Figures/af_figures/pfs4845_AF1.png, width=800, height=800, dpi=800, ray=1

rotate y, 180

png ../Pymol Figures/af_figures/pfs4845_AF2.png, width=800, height=800, dpi=800, ray=1