## Finalize wattersons pictures

# background white and alphachecker off
bg_color white
set ray_opaque_background, off

load ../tajima_d_models/ama1_full_norsa_trim_pi_adj.pdb
load ../tajima_d_models/help_scale5.pdb

# set view (orient + zoom center, 60)
set_view (\
    -0.600028574,   -0.126176745,   -0.789965272,\
     0.079115048,    0.973282635,   -0.215549886,\
     0.796056867,   -0.191834226,   -0.574014902,\
     0.000000000,    0.000000000, -352.726074219,\
     0.839885712,   -1.681072235,   -1.454235077,\
   280.726074219,  424.726074219,  -20.000000000 )

# run quick_color2
as surface

# turn off help_scale2
hide everything, help_scale5

#Select all residues with a mapped data value. Can change the default 'no-value'
#option when writing to pdb b factor using biostructmap if needed.
select nonzeros, b < 0 
select tz, not (b > 0 | b < 0)

color grey
#color white, tz

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

png ../Pymol Figures/ama1_pi1_full.png, width=800, height=800, dpi=800, ray=1

rotate y, 180

png ../Pymol Figures/ama1_pi2_full.png, width=800, height=800, dpi=800, ray=1