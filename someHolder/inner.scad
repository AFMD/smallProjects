//Outer block
height=50;
width=50;
depth=19;

//substrate slot

sHeight=12.5;
sWidth=12.5;
sDepth=9;

//Ridge
ridDepth=0.074*25.4; //As per workshop recommendation
ridInner=25.1/2; // radius for o-ring #120
ridCross=2.6; //Cross section

//canister
rRadius=ridInner+ridCross+2;
rDepth=sDepth+2;

//accuracy
n=250;




//translate([rDepth/2+ridDepth/2,0,0])rotate([0,90,0])cylinder(ridDepth,rRadius-2,rRadius-2,center=true);}


   
   
   //Substrate slot in canister 
    //Inner Ridge -Front
   dis=50;
   //Substrate slot in canister 
    translate([ridDepth/4,dis,0])rotate([0,90,0])difference(){cube([2*rRadius,2*rRadius,rDepth-ridDepth/2],center=true);
    rotate([0,90,0])cube([sDepth,sWidth,sHeight],center=true);
    rotate([0,90,0])translate([-sDepth/2,0,0])cube([2,sHeight-1.5,sWidth-1.5],center=true);
    rotate([0,90,0])translate([sDepth/2,0,0])cube([2,sHeight,sWidth],center=true);}
    
    difference(){translate([-rDepth/2+ridDepth/2,dis,0])rotate([0,90,0])cube([2*rRadius,2*rRadius,ridDepth],center=true);
       translate([-rDepth/2+ridDepth/2,dis,0])rotate([0,90,0])cylinder(ridDepth,ridInner+ridCross,ridInner+ridCross,center=true,$fn=n);}
    //inner ridge front
      difference(){translate([-rDepth/2+ridDepth/2,dis,0])rotate([0,90,0])cylinder(ridDepth,ridInner,ridInner,center=true,$fn=n);;
       translate([-rDepth/2+ridDepth/2,dis,0])rotate([0,0,0])cube([ridDepth+1,sHeight,sWidth],center=true);}
       


    
    
 
