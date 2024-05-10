## Finalize wattersons pictures

# background white and alphachecker off
bg_color white
set ray_opaque_background, off

load ../tajima_d_models/msp119_norsa_trim_tajima.pdb
load ../tajima_d_models/help_scale.pdb

# set view (orient + zoom center, 60)
set_view (\
    -0.564158559,    0.785217762,    0.255261004,\
     0.331755161,   -0.067528039,    0.940945625,\
     0.756084442,    0.615526617,   -0.222403511,\
     0.000000000,    0.000000000, -340.276916504,\
    -2.851646423,   -1.809741974,    0.455833435,\
   268.276916504,  412.276916504,  -20.000000000 )

# run quick_color2
as surface

# turn off help_scale
hide everything, help_scale

#Select all residues with a mapped data value. Can change the default 'no-value'
#option when writing to pdb b factor using biostructmap if needed.
select nonzeros, b < 0 | b > 0

color grey

spectrum b, blue_white_magenta, selection="nonzeros"

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

png ../Pymol Figures/msp119_taj1.png, width=800, height=800, dpi=800, ray=1

rotate y, 180

png ../Pymol Figures/msp119_taj2.png, width=800, height=800, dpi=800, ray=1