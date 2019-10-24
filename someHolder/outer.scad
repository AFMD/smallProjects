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

//Outer holder
difference(){translate([ridDepth/2,0,0])cube([rDepth+ridDepth,width,height],center=true);
    rotate([0,90,0]){cube([2*rRadius+1.5,2*rRadius+1.5,rDepth],center=true);;
    }
    //Outer ridge - back
    translate([rDepth/2+ridDepth/2,0,0])rotate([0,90,0])cylinder(ridDepth,ridInner+ridCross,ridInner+ridCross,center=true,$fn=n);
    
    //screw holes
    translate([-rDepth/2,0.4*width,0.4*height])rotate([0,90,0])cylinder(3,2,2,center=true,$fn=n);
     translate([-rDepth/2,-0.4*width,0.4*height])rotate([0,90,0])cylinder(3,2,2,center=true,$fn=n);
     translate([-rDepth/2,0.4*width,-0.4*height])rotate([0,90,0])cylinder(3,2,2,center=true,$fn=n);
     translate([-rDepth/2,-0.4*width,-0.4*height])rotate([0,90,0])cylinder(3,2,2,center=true,$fn=n);
    }
    //back stop
difference(){translate([ridDepth/2+rDepth/2+1.5,0,0])cube([3,width,height],center=true);
    translate([ridDepth/2+rDepth/2+1.5,0,0])cube([4,sHeight-1.5,sWidth-1.5],center=true);
}
    //Inner Ridge -Front
difference(){translate([rDepth/2+ridDepth/2,0,0])rotate([0,90,0])cylinder(ridDepth,ridInner,ridInner,center=true,$fn=n);
    translate([rDepth/2+ridDepth/2,0,0])cube([ridDepth+1,sHeight-1.5,sWidth-1.5],center=true);}
