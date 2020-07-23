include <meta.scad>;

$fn=50;
deviceDim=28; //mm, x,y dimension of substrate, pcb
//deviceDim=30; //mm, x,y dimension of substrate, pcb
//deviceDim=25; //mm, x,y dimension of substrate, pcb, DBG
capFit=0.4; //mm, trace a circle of this radius around the crown pieces to define the cap cutouts
//capX = 39;//mm
//capY = 60;//mm

vDeviceDim=deviceDim+deviceFudge;
xyLen=vDeviceDim+wallT*2;

capX = xyLen;
capY = xyLen;

windowDim = deviceDim-4*2;
capT = 5;
rotate([180,0,0])
difference(){
    linear_extrude(height=capT){
        difference(){
            square([capX,capY], center=true);
            minkowski(){
              projection(cut = true) translate([0,0,-totalHeight/2+1]) base(deviceDim);
              circle(capFit);
            }
            square([windowDim,windowDim], center=true);
        }
    }
    
    linear_extrude(height=capT,scale=windowDim/deviceDim){
        square([deviceDim,deviceDim], center=true);
    }
}