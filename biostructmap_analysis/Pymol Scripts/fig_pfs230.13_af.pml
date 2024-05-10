## show AF confidence

# background white and alphachecker off
bg_color white
set ray_opaque_background, off

load ../biostructmap_input/pfs230_seg13.pdb
load ../tajima_d_models/help_scale_AF.pdb

# set view (orient + zoom center, 60)
set_view (\
     0.685389698,    0.412666976,    0.599955797,\
    -0.449356616,    0.888020277,   -0.097460903,\
    -0.572991848,   -0.202795416,    0.794074595,\
     0.000000000,    0.000000000, -340.276916504,\
    -0.956398010,    0.466567993,   -0.133190155,\
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

png ../Pymol Figures/af_figures/pfs230_seg13_af1.png, width=800, height=800, dpi=800, ray=1

rotate y, 180

png ../Pymol Figures/af_figures/pfs230_seg13_af2.png, width=800, height=800, dpi=800, ray=1