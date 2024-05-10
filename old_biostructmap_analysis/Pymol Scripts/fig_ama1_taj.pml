## Finalize wattersons pictures

# background white and alphachecker off
bg_color white
set ray_opaque_background, off

load ../tajima_d_models/ama1_trim_tajima.pdb
load ../tajima_d_models/help_scale.pdb

# set view (orient + zoom center, 60)
set_view (\
    -0.880948722,   -0.116181947,    0.458725601,\
    -0.303812057,    0.882078469,   -0.360046148,\
    -0.362802237,   -0.456547588,   -0.812367082,\
    -0.000000000,    0.000000000, -340.276916504,\
    12.691177368,  -31.270515442,   36.808609009,\
   268.276916504,  412.276916504,  -20.000000000 )

# run quick_color2
as surface

# turn off help_scale2
hide everything, help_scale

#Select all residues with a mapped data value. Can change the default 'no-value'
#option when writing to pdb b factor using biostructmap if needed.
select nonzeros, b < 0 | b > 0

color white

spectrum b, selection=nonzeros

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

png ../Pymol Figures/ama1_taj1.png, width=800, height=800, dpi=800, ray=1

rotate y, 180

png ../Pymol Figures/ama1_taj2.png, width=800, height=800, dpi=800, ray=1