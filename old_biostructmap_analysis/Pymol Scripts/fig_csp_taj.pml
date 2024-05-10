## Finalize wattersons pictures

# background white and alphachecker off
bg_color white
set ray_opaque_background, off

load ../tajima_d_models/csp_trim_tajima.pdb
load ../tajima_d_models/help_scale.pdb

# set view (orient + zoom center, 60)
set_view (\
     0.359280556,   -0.129196391,    0.924243331,\
    -0.411691338,    0.866851330,    0.281210303,\
    -0.837512851,   -0.481536299,    0.258253783,\
     0.000000000,    0.000000000, -340.276916504,\
    32.065567017,   63.320732117,  -88.451034546,\
   268.276916504,  412.276916504,  -20.000000000 )

# run quick_color
as surface

# turn off help_scale
hide everything, help_scale

#Select all residues with a mapped data value. Can change the default 'no-value'
#option when writing to pdb b factor using biostructmap if needed.
select nonzeros, b < 0 | b > 0

color white

spectrum b, selection="nonzeros"

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

png ../Pymol Figures/csp_taj1.png, width=800, height=800, dpi=800, ray=1

rotate y, 180

png ../Pymol Figures/csp_taj2.png, width=800, height=800, dpi=800, ray=1



