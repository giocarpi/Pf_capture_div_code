## show AF confidence

# background white and alphachecker off
bg_color white
set ray_opaque_background, off

load ../biostructmap_input/pfs47_AF_cleaned_short.pdb
load ../tajima_d_models/help_scale_AF.pdb

# set view (orient + zoom center, 60)
set_view (\
     0.865265369,   -0.152295053,    0.477621228,\
    -0.372458190,    0.442402273,    0.815815628,\
    -0.335545391,   -0.883790910,    0.326071829,\
     0.000000000,    0.000000000, -340.276916504,\
    -1.416248322,   -0.204406738,    0.981872559,\
   268.276916504,  412.276916504,  -20.000000000 )

# run quick_color2
as cartoon

# turn off help_scale2
hide everything, help_scale_AF

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

png ../Pymol Figures/af_figures/pfs47_short_af1.png, width=800, height=800, dpi=800, ray=1

rotate y, 180

png ../Pymol Figures/af_figures/pfs47_short_af2.png, width=800, height=800, dpi=800, ray=1