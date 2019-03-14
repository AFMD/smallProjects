// design for a sphere holder for Niki's quartz
// grey@mutovis.com
// 12 March 2019

include <../../mutovis/logo/logo_modules.scad>
$fn=100;


// all dims in mm
major_cap = [40, 50, 3]; // holder outer dims
major_base = [50, 40, 3]; // holder outer dims
substrate = [15, 15, 1]; // substrate dimensions
substrate_fudge = 0.6; // manufacturing fudge factor
shelf_thickness = 0.5;
shelf_width = 2;
window = [substrate[0]-shelf_width, substrate[1]-shelf_width];
substrate_fudged = [substrate[0]+substrate_fudge, substrate[1]+substrate_fudge];

module base(outer) {
    difference() {
        translate([-outer[0]/2, -outer[1]/2, 0]) cube(outer);
        translate([-window[0]/2, -window[1]/2, -1]) cube([window[0], window[0],outer[2]*2]);
        translate([-substrate_fudged[0]/2, -substrate_fudged[1]/2, shelf_thickness]) cube([substrate_fudged[0], substrate_fudged[1], outer[2]*2]);
        }
}

module cap(outer) {
    difference() {
        union(){
            translate([-outer[0]/2, -outer[1]/2, 0]) cube(outer);
            translate([-substrate[0]/2, -substrate[1]/2, 0]) cube([substrate[0], substrate[1], 2*outer[2]-shelf_thickness-substrate[2]]);
        }
        translate([-window[0]/2, -window[1]/2, -1]) cube([window[0], window[0],outer[2]*2]);
    }
}

module fancy(){
    logo_scale = 0.7;
    text_scale = 0.456;
    mid_offset = (major_cap[1]/2 - substrate[1]/2)/2+substrate[1]/2;
    linear_extrude(height = major_cap[2]){
        translate([0,mid_offset,0]) scale([-logo_scale, logo_scale, 1]) logo();
        translate([0,-mid_offset,0]) scale([-text_scale, text_scale, 1]) logo_text();
    }
}

// STL1:
//base(major_base);

// STL2:
difference(){
    cap(major_cap);
    fancy();
}

// STL3:
//fancy();