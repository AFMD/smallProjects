// design for a sphere holder for Niki's quartz
// grey@mutovis.com
// 12 March 2019

// all dims in mm
major = [40, 40, 3]; // holder outer dims
substrate = [15, 15, 1]; // substrate dimensions
substrate_fudge = 0.2; // manufacturing fudge factor
shelf_thickness = 0.5;
shelf_width = 2;
window = [substrate[0]-shelf_width, substrate[1]-shelf_width];
substrate_fudged = [substrate[0]+substrate_fudge, substrate[1]+substrate_fudge];

module base() {
    difference() {
        translate([-major[0]/2, -major[1]/2, 0]) cube(major);
        translate([-window[0]/2, -window[1]/2, -1]) cube([window[0], window[0],major[2]*2]);
        translate([-substrate_fudged[0]/2, -substrate_fudged[1]/2, shelf_thickness]) cube([substrate_fudged[0], substrate_fudged[1], major[2]*2]);
        }
}

module cap() {
    difference() {
        union(){
            translate([-major[0]/2, -major[1]/2, 0]) cube(major);
            translate([-substrate[0]/2, -substrate[1]/2, 0]) cube([substrate[0], substrate[1], 2*major[2]-shelf_thickness-substrate[2]]);
        }
        translate([-window[0]/2, -window[1]/2, -1]) cube([window[0], window[0],major[2]*2]);
        
    }
}

base();
//cap();