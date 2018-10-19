glass_thickness = 1.25 + 0.1;
removal_notch = 5;
glass_dim = 28;
glass_fudge = 0.2;
pocket_dim = glass_dim + glass_fudge;

//$fn = 100;
extraXY = 50;

x = pocket_dim + extraXY;
y = pocket_dim + extraXY;
z = 5;

slot_x_position = extraXY/2;

slot_depth = z - glass_thickness;

slot_width = 4;



difference(){
    cube([x,y,z]);
    
    // first slot
    translate([slot_x_position,0,glass_thickness]) cube([slot_width,y,slot_depth]);
    
    // middle slot
    translate([x/2-slot_width/2,0,glass_thickness]) cube([slot_width,y,slot_depth]);
    
    // third slot
    translate([x-slot_x_position-slot_width,0,glass_thickness]) cube([slot_width,y,slot_depth]);
    
    // glass pocket
    translate([(x-pocket_dim)/2,(y-pocket_dim)/2,0]) cube([pocket_dim,pocket_dim,glass_thickness]);
    
    // glass removal notch
    
    translate([(x-pocket_dim)/2,(extraXY+pocket_dim)/2,0]) cylinder(glass_thickness,removal_notch,removal_notch);
    
}
