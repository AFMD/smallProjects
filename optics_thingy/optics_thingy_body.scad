// units in mm
// optics light tube for Akash
$fn = 40;

// editable
post_diameter_nominal = 12.7;
post_d_fudge = -0.4;

post_height = 90;

platexy = 50;
plate_t = 5;
inner_d = 10;


colm_diameter_nominal = 12.7;
colm_d_fudge = 0.3;
colm_h = 15;

// collimator support cylinder
support_cyl_d_offset = 5;
support_cyl_h_offset = 3;

mount_hole_d = 6.5;
hole_offset = 15;

// calculated
post_d = post_diameter_nominal + post_d_fudge;
colm_d = colm_diameter_nominal + colm_d_fudge;
holexy = platexy - hole_offset;
support_cyl_d = colm_h + support_cyl_d_offset;
support_cyl_h = colm_h + support_cyl_h_offset;

module body(){
    // geometry
    difference(){
        union(){
            cylinder(d=post_d, h=post_height);
            cylinder(d=support_cyl_d, h=support_cyl_h);
            translate([0,0,plate_t/2]) cube([platexy, platexy, plate_t],center=true);
        }
        cylinder(d=inner_d, h=post_height);
        cylinder(d=colm_d, h=colm_h);
        
        translate([-holexy/2,-holexy/2,0]) cylinder(d=mount_hole_d, h=plate_t);
        translate([ holexy/2,-holexy/2,0]) cylinder(d=mount_hole_d, h=plate_t);
        translate([-holexy/2, holexy/2,0]) cylinder(d=mount_hole_d, h=plate_t);
        translate([ holexy/2, holexy/2,0]) cylinder(d=mount_hole_d, h=plate_t);
    }
}

body();