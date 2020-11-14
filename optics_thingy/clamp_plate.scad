include <optics_thingy_body.scad>

fiber_clearance_d = 10;

// geometry

translate([0,0,-plate_t]) difference(){
    translate([0,0,plate_t/2]) cube([platexy, platexy, plate_t],center=true);
    cylinder(d=fiber_clearance_d, h=colm_h);
    
    translate([-holexy/2,-holexy/2,0]) cylinder(d=mount_hole_d, h=plate_t);
    translate([ holexy/2,-holexy/2,0]) cylinder(d=mount_hole_d, h=plate_t);
    translate([-holexy/2, holexy/2,0]) cylinder(d=mount_hole_d, h=plate_t);
    translate([ holexy/2, holexy/2,0]) cylinder(d=mount_hole_d, h=plate_t);
}

//body();