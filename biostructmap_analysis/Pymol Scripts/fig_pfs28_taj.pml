## Finalize wattersons pictures

# background white and alphachecker off
bg_color white
set ray_opaque_background, off

load ../tajima_d_models/pfs28_norsa_trim_tajima.pdb
load ../tajima_d_models/help_scale.pdb

# set view (orient + zoom center, 60)
set_view (\
     0.275250256,    0.935057282,    0.223394483,\
    -0.925202012,    0.194503546,    0.325836748,\
     0.261225015,   -0.296371669,    0.918654084,\
     0.000000000,    0.000000000, -340.276916504,\
    14.997909546,   -1.630996704,    7.737716675,\
   268.276916504,  412.276916504,  -20.000000000 )

# run quick_color2
as surface

# turn off help_scale2
hide everything, help_scale

#Select all residues with a mapped data value. Can change the default 'no-value'
#option when writing to pdb b factor using biostructmap if needed.

select nonzeros, b < 0 | b > 0

color grey

spectrum b, blue_white_magenta, selection=nonzeros

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

png ../Pymol Figures/pfs28_taj1.png, width=800, height=800, dpi=800, ray=1

rotate y, 180

png ../Pymol Figures/pfs28_taj2.png, width=800, height=800, dpi=800, ray=1