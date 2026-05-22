from math import cos, sin, radians
import cadquery as cq
from cadquery.vis import show
#from cadquery.plugins import Xutil
#from cadquery.plugins import shapex

h1= 2
r1, r1a =30 ,6
r2= 23.1
r3=3.5
dx= r3*cos(radians(30))
dy= r3*sin(radians(30))
h2= 90
r4, r5= 26, 28

# helper function
def polar_1array(obj, num_copies, angle_deg):
    result = cq.Workplane("XY")
    for i in range(num_copies):
        result = result.union(obj.rotate((0, 0, 0), (0, 0, 1), angle_deg * i))
    return result

def Sector(radius, angle, thickness):
    """
    Creates a sector (pie slice) shape.
    Parameters:
    - radius: float - radius of the sector
    - angle: float - angle of the sector in degrees
    - thickness: float - thickness of the extruded sector
    Returns:
    - CadQuery workplane object containing the sector
    """

    # Convert angle to radians (half angle for the calculation)
    alpha = radians(angle/2)

    # Calculate points for the sector
    a = radius * sin(alpha)
    b = radius * cos(alpha)

    # Create the sector shape
    rs = (cq.Workplane("XY")

    .moveTo(0, 0) # Start at origin
    .lineTo(a, b) # Line to first point
    .threePointArc((0, radius), (-a, b)) # Arc through top point
    .lineTo(0, 0) # Line back to origin
    .close() # Close the sketch
    .extrude(thickness) # Extrude to create 3D shape
    )

    show(rs)
    return rs

def annular(r1,r2,angle,t1):

    #r2> r1
    bound1= Sector(r1,angle,t1)
    bound2= Sector(r2,angle,t1)
    rs= bound2.cut(bound1)
    return rs

#wp= cq.Workplane("XZ")
wp1= cq.Workplane("XY")
bas1= wp1.circle(r1).extrude(h1).translate((0,0,0))
bas2= wp1.circle(r1a).extrude(h1)
sec= Sector(r2,120,h1).rotate((0,0,0), (0,0,1), 0)
sec= sec.translate((0,dy,0))
g_sec= polar_1array(sec, 3, 120)
bind1= bas1.cut(g_sec)
bind1= bind1.union(bas2)
bind2= bind1.translate((0,0,h2))
wing= annular(r4,r5,26,h2)
g_wing= polar_1array(wing, r4, 90)
frame= bind1.union(bind2).union(g_wing)

#show
show(frame)
show(g_wing)

