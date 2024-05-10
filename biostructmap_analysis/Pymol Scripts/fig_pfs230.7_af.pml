## Finalize wattersons pictures

# background white and alphachecker off
bg_color white
set ray_opaque_background, off

load ../biostructmap_input/pfs230_seg7.pdb
load ../tajima_d_models/help_scale_AF.pdb

# set view (orient + zoom center, 60)
set_view (\
     0.611450315,    0.683827519,   -0.398131102,\
    -0.319515944,   -0.246931612,   -0.914841294,\
    -0.723904729,    0.686589181,    0.067507371,\
     0.000000000,    0.000000000, -340.276916504,\
     1.095432281,   -3.112777710,   -3.050315857,\
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

png ../Pymol Figures/af_figures/pfs230_seg7_af1.png, width=800, height=800, dpi=800, ray=1

rotate y, 180

png ../Pymol Figures/af_figures/pfs230_seg7_af2.png, width=800, height=800, dpi=800, ray=1