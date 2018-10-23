glass_thickness = 1.1;
glass_z_fudge = 0.1;
pocket_thickness = glass_thickness + glass_z_fudge;
active_area_buffer_thickness = 0.25;

glass_dim = 38;
glass_fudge = 0.5;
pocket_dim = glass_dim + glass_fudge;
ridge_dim = 10;

extraXY = 50;

x = pocket_dim + extraXY;
y = pocket_dim + extraXY;
z = 3;

slot_depth = z - glass_thickness;
slot_width = 4;
wide_slot_width = 25;
blocker_width = (pocket_dim - (3*slot_width))/4;
slot_x_position = extraXY/2;
removal_notch = blocker_width;
// Ridge for structural stability
difference(){
    cube([x,y,z]);
    
    translate([ridge_dim,ridge_dim,0])
    cube([x-2*ridge_dim,y-2*ridge_dim,z]);
}
difference(){
    cube([x,y,z]);
    
    // first slot
    translate([extraXY/2 - wide_slot_width+slot_width,0,glass_thickness]) cube([wide_slot_width,y,slot_depth]);
    
    // middle slot
    translate([x/2-slot_width/2,0,glass_thickness]) cube([slot_width,y,slot_depth]);
    
    // third slot
    translate([x-slot_x_position-slot_width,0,glass_thickness]) cube([wide_slot_width,y,slot_depth]);
    
    // glass pocket
    translate([(x-pocket_dim)/2,(y-pocket_dim)/2,0]) cube([pocket_dim,pocket_dim,pocket_thickness]);
    
    // active area buffer
    translate([(x-pocket_dim)/2+slot_width/2,(y-pocket_dim)/2+slot_width/2,0]) cube([pocket_dim-slot_width,pocket_dim-slot_width,glass_thickness+active_area_buffer_thickness]);
    
    // glass removal notch
    
    translate([x/2+(slot_width/2)+blocker_width,extraXY/2+pocket_dim,0]) cylinder(glass_thickness,removal_notch,removal_notch);
    
}
