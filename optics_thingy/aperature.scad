include <optics_thingy_body.scad>

aper_t = 3;
aperxy = 5;

translate([0,0,post_height]) difference(){
    cylinder(d=post_d,h=aper_t);
    translate([0,0,aper_t/2]) cube([aperxy, aperxy, aper_t], center=true);
}
//body();