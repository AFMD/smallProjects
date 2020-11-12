// units in mm
// optics light tube for Akash

// editable
post_diameter_nominal = 25;
post_d_fudge = -0.1;

post_height = 60;

platexy = 50;
plate_t = 10;
inner_d = 10;

colm_diameter_nominal = 12;
colm_d_fudge = 0.1;
colm_h = 12.7;

mount_hole_d = 5;
hole_offset = 15;



// calculated
post_d = post_diameter_nominal + post_d_fudge;
colm_d = colm_diameter_nominal + colm_d_fudge;
holexy = platexy - hole_offset;

// geometry
difference(){
    union(){
        cylinder(d=post_d, h=post_height);
        translate([0,0,plate_t/2]) cube([platexy, platexy, plate_t],center=true);
    }
    cylinder(d=inner_d, h=post_height);
    cylinder(d=colm_d, h=colm_h);
    
    translate([-holexy/2,-holexy/2,0]) cylinder(d=mount_hole_d, h=plate_t);
    translate([ holexy/2,-holexy/2,0]) cylinder(d=mount_hole_d, h=plate_t);
    translate([-holexy/2, holexy/2,0]) cylinder(d=mount_hole_d, h=plate_t);
    translate([ holexy/2, holexy/2,0]) cylinder(d=mount_hole_d, h=plate_t);

}