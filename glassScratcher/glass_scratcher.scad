glass_thickness = 1.1;
no_scratch_z = 0.25;
pocket_thickness = glass_thickness + no_scratch_z;

glass_dim = 30;
glass_xy_fudge = 0.5;
pocket_xy = glass_dim + glass_xy_fudge;

// for mechanical robustness
outer_rim_width = 15;

extraXY = 50;

// main cube outer dims
x = pocket_xy + extraXY;
y = pocket_xy + extraXY;
z = 4; // for mechanical robustness

slot_width = 1.5;
side_opening_width = 5.75;

side_pocket_x = y/2-outer_rim_width-glass_dim/2+side_opening_width;

difference(){
    cube([x,y,z]);
    
    // first slot
    translate([outer_rim_width, outer_rim_width, glass_thickness]) cube([side_pocket_x, y-2*outer_rim_width, z-glass_thickness]);
    
    // middle slot
    translate([x/2-slot_width/2,outer_rim_width, glass_thickness]) cube([slot_width,y-2*outer_rim_width, z-glass_thickness]);
    
    // third slot
    translate([x-outer_rim_width-side_pocket_x, outer_rim_width, glass_thickness]) cube([side_pocket_x, y-2*outer_rim_width, z-glass_thickness]);
    
    // glass pocket
    translate([(x-pocket_xy)/2,(y-pocket_xy)/2,0]) cube([pocket_xy, pocket_xy, pocket_thickness]);
}