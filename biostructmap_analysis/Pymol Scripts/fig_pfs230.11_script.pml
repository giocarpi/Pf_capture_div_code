## Finalize wattersons pictures

# background white and alphachecker off
bg_color white
set ray_opaque_background, off

load ../tajima_d_models/pfs230.11_norsa_trim_pi_adj.pdb
load ../tajima_d_models/help_scale5.pdb

# set view (orient + zoom center, 60)
set_view (\
     0.780339181,   -0.563463211,    0.271256775,\
    -0.305936754,    0.034337770,    0.951432407,\
    -0.545411468,   -0.825427294,   -0.145588994,\
     0.000000000,    0.000000000, -340.276916504,\
     4.739570618,   -3.164241791,   -2.903465271,\
   268.276916504,  412.276916504,  -20.000000000 )

# run quick_color2
as surface

# turn off help_scale
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

png ../Pymol Figures/pfs230_seg11_pi1.png, width=800, height=800, dpi=800, ray=1

rotate y, 180

png ../Pymol Figures/pfs230_seg11_pi2.png, width=800, height=800, dpi=800, ray=1