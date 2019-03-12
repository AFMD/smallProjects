// design for a sphere holder for Niki's quartz
// grey@mutovis.com
// 12 March 2019

// all dims in mm
major_cap = [40, 50, 3]; // holder outer dims
major_base = [50, 40, 3]; // holder outer dims
substrate = [15, 15, 1]; // substrate dimensions
substrate_fudge = 0.4; // manufacturing fudge factor
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

base(major_base);
//cap(major_cap);