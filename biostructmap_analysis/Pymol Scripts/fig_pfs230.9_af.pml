## Finalize wattersons pictures

# background white and alphachecker off
bg_color white
set ray_opaque_background, off

load ../biostructmap_input/pfs230_seg9.pdb
load ../tajima_d_models/help_scale_AF.pdb

# set view (orient + zoom center, 60)
set_view (\
     0.732946098,   -0.050562337,    0.678405106,\
    -0.222988755,    0.924282193,    0.309803993,\
    -0.642702162,   -0.378346324,    0.666174114,\
     0.000000000,    0.000000000, -340.276916504,\
    -3.020706177,   -2.551059723,   -0.177989960,\
   268.276916504,  412.276916504,  -20.000000000 )

# run quick_color2
as cartoon

# turn off help_scale
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

png ../Pymol Figures/af_figures/pfs230_seg9_af1.png, width=800, height=800, dpi=800, ray=1

rotate y, 180

png ../Pymol Figures/af_figures/pfs230_seg9_af2.png, width=800, height=800, dpi=800, ray=1