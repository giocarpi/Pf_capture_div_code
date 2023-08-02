## Finalize wattersons pictures

# background white and alphachecker off
bg_color white
set ray_opaque_background, off

load ../tajima_d_models/rh5_trim_pi_adj.pdb
load ../tajima_d_models/help_scale5.pdb

# set view (orient + zoom center, 60)
set_view (\
     0.281602055,   -0.871823132,    0.400780141,\
    -0.452443331,    0.247689858,    0.856705904,\
    -0.846165121,   -0.422580391,   -0.324700534,\
     0.000000000,    0.000000000, -340.276916504,\
    27.980653763,   49.052391052,   75.772521973,\
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

png ../Pymol Figures/rh5_pi1.png, width=800, height=800, dpi=800, ray=1

rotate y, 180

png ../Pymol Figures/rh5_pi2.png, width=800, height=800, dpi=800, ray=1