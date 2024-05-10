## show AF confidence

# background white and alphachecker off
bg_color white
set ray_opaque_background, off

load ../biostructmap_input/pfs230_seg3.pdb
load ../tajima_d_models/help_scale_AF.pdb

# set view (orient + zoom center, 60)
set_view (\
     0.373656392,    0.918478727,   -0.129529059,\
     0.131368607,   -0.190637589,   -0.972830713,\
    -0.918217421,    0.346488357,   -0.191892222,\
     0.000000000,    0.000000000, -340.276916504,\
     0.172626495,    2.064426422,   -4.639022827,\
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

png ../Pymol Figures/af_figures/pfs230_seg3_af1.png, width=800, height=800, dpi=800, ray=1

rotate y, 180

png ../Pymol Figures/af_figures/pfs230_seg3_af2.png, width=800, height=800, dpi=800, ray=1