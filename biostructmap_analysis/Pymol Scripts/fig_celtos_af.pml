## show AF confidence

# background white and alphachecker off
bg_color white
set ray_opaque_background, off

load ../biostructmap_input/celtos_AF_cleaned.pdb
load ../tajima_d_models/help_scale_AF.pdb

# set view (orient + zoom center, 60)
set_view (\
    -0.555667996,    0.425812721,   -0.714084566,\
     0.145033970,    0.895365655,    0.421052814,\
     0.818656385,    0.130399033,   -0.559283316,\
     0.000000000,    0.000000000, -340.276916504,\
    -7.038619995,    4.020566940,    1.632442474,\
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

png ../Pymol Figures/af_figures/celtos_af1.png, width=800, height=800, dpi=800, ray=1

rotate y, 180

png ../Pymol Figures/af_figures/celtos_af2.png, width=800, height=800, dpi=800, ray=1