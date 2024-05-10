## Finalize wattersons pictures

# background white and alphachecker off
bg_color white
set ray_opaque_background, off

load ../tajima_d_models/pfs230_seg1_trim_pi_adj.pdb
load ../tajima_d_models/help_scale5.pdb

# set view (orient + zoom center, 60)
set_view (\
     0.527006626,    0.088863380,    0.845202506,\
    -0.086744592,    0.994948804,   -0.050519895,\
    -0.845422626,   -0.046692427,    0.532052994,\
     0.000000000,    0.000000000, -340.276916504,\
    -0.525707245,    2.421504974,    2.701141357,\
   268.276916504,  412.276916504,  -20.000000000 )

# run quick_color2
as surface

# turn off help_scale2
hide everything, help_scale5

#Select all residues with a mapped data value. Can change the default 'no-value'
#option when writing to pdb b factor using biostructmap if needed.
select nonzeros, b < 0 | b > 0

color white

spectrum b, blue_red, selection="nonzeros"

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

png ../Pymol Figures/pfs230_seg1_pi1.png, width=800, height=800, dpi=800, ray=1

rotate y, 180

png ../Pymol Figures/pfs230_seg1_pi2.png, width=800, height=800, dpi=800, ray=1