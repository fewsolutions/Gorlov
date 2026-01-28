def generator(af1,                      #float, max camber of airfoil (first number of NACA 4 digit code)
              af2,                      #float, location of max camber of airfoil (second number of NACA code)
              af34,                     #max thickness (last two numbers of NACA code)  CANNOT BE 0
              numBlades,                #int, number of blades (max 9 for centreRadius=4) 
              airfoilLength,            #float, length of the airfoil cross-section
              airfoilRotation,          #float, rotation angle of airfoil in degrees in relation to tangent of circle 
    ):
    
    diameter=30              #float, turbine outer diameter 
    length=40                #float, turbine total length/height 
    outerRingThickness=3     #float, thickness of outer ring on end caps 
    spokesThickness=3        #float, thickness of spokes on end caps 
    centreRadius=3           #float, radius of center hub 
    endsThickness=2          #float, thickness of end caps 
    airfoilNoPoints=200      #int, resolution of NACA airfoil cross-section
    filename="turbine"       #str, export filename

    if numBlades > 9:
        raise ValueError("Too many blades for this centreRadius")

    import math
    import cadquery as cq
    from naca import naca

    radius = diameter/2
    bladesradius = radius - outerRingThickness/2
    middleCircleRadius = radius - outerRingThickness
    centre = (0,0)
    

    #--- get inner circle points ---
    def getCoordinate(pointNo, direction, numBlades, radius, spokesThickness):
        theta = ((2 * math.pi) / numBlades) * pointNo  # angle in radians
        offset = spokesThickness / 2
        
        t = math.sqrt(radius**2 - offset**2)

        # Move along tangent by offset either clockwise or counterclockwise
        if direction == "ccw":
            x = centre[0] + t * math.cos(theta) - offset * math.sin(theta)
            y = centre[1] + t * math.sin(theta) + offset * math.cos(theta)
        elif direction == "cw":
            x = centre[0] + t * math.cos(theta) + offset * math.sin(theta)
            y = centre[1] + t * math.sin(theta) - offset * math.cos(theta)
        elif direction == "mid":
            x = centre[0] + radius * math.cos(theta)  #original point (x)
            y = centre[1] + radius * math.sin(theta)  #original point (y)
        return (round(x, 15), round(y, 15))

    #--- create end caps ---
    sketch = (cq.Sketch("XY"))
    for i in range(numBlades):
        sketch = sketch.segment(getCoordinate(i, "ccw", numBlades, centreRadius, spokesThickness), getCoordinate(i, "ccw", numBlades, middleCircleRadius, spokesThickness)).segment(getCoordinate(i, "cw", numBlades, centreRadius, spokesThickness), getCoordinate(i, "cw", numBlades, middleCircleRadius, spokesThickness)).arc(getCoordinate(i, "ccw", numBlades, middleCircleRadius, spokesThickness), getCoordinate(i*2+1, "mid", numBlades*2, middleCircleRadius, spokesThickness), getCoordinate(i + 1, "cw", numBlades, middleCircleRadius, spokesThickness)).arc(getCoordinate(i, "ccw", numBlades, centreRadius, spokesThickness), getCoordinate(i*2+1, "mid", numBlades*2, centreRadius, spokesThickness), getCoordinate(i + 1, "cw", numBlades, centreRadius, spokesThickness))
    # draw outside circle workaround using two half circles
    sketch = sketch.arc(getCoordinate(0, "mid", 4, radius, spokesThickness), getCoordinate(0*2+1, "mid", 8, radius, spokesThickness), getCoordinate(2, "mid", 4, radius, spokesThickness)).arc(getCoordinate(0, "mid", 4, radius, spokesThickness), getCoordinate(0*2-1, "mid", 8, radius, spokesThickness), getCoordinate(2, "mid", 4, radius, spokesThickness)).assemble()

    # extrude sketch down to make 3D end cap
    model = cq.Workplane("XY").placeSketch(sketch).extrude(-endsThickness)

    #--- get airfoil points and scale correctly ---
    airfoilPoints = naca(af1, af2, af34, n=airfoilNoPoints)
    airfoilPoints = [(x * airfoilLength, y * airfoilLength) for x, y in airfoilPoints]

    #--- build helix ---
    pitch = (length - endsThickness * 2) * numBlades
    height = length - endsThickness * 2
    wire = cq.Wire.makeHelix(pitch=pitch, height=height, radius=bladesradius)
    helix = cq.Workplane(obj=wire)
    for i in range(numBlades):
        x, y = getCoordinate(i, "mid", numBlades, bladesradius, spokesThickness)
        # Angle of the radius to the placement point
        theta = ((2 * math.pi) / numBlades) * i
        # Tangent direction (CCW) is radius angle + 90 degrees; include user airfoilRotation (degrees)
        tangent_angle = theta + math.pi / 2 + math.radians(airfoilRotation)
        ca, sa = math.cos(tangent_angle), math.sin(tangent_angle)
        # Rotate the scaled airfoil points around origin to align with tangent
        rotated_airfoil = [(px * ca - py * sa, px * sa + py * ca) for (px, py) in airfoilPoints]
        # Build blade profile at absolute (x,y) and sweep along helix
        blade = (
            cq.Workplane('XY')
            .center(x, y)
            .spline(rotated_airfoil)
            .close()
            .sweep(helix, isFrenet=True)
        )
        model = model.add(blade)

    # add top cap to model
    topCap = (
        cq.Workplane("XY")
        .workplane(offset=height)
        .placeSketch(sketch)
        .extrude(endsThickness)
        )
    model = model.add(topCap)

    # Export model as .stl
    #cq.exporters.export(model, f"{filename}.stl")
    # Export model as .step
    cq.exporters.export(model, f"{filename}.step")