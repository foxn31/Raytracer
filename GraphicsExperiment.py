# GraphicsExperiment.py
# Nick Fox
# Experimental program playing with graphics, the goal being to make a raytracer
# Uses graphics.py by John Zelle

`from graphics import *
`
def main():
    win = GraphWin("My Circle", 1000, 800)
    c = Rectangle(Point(50,50),Point(100,100))
    c.draw(win)

    tri = Tri(((10,20),(20,20),(10,10)))
    point = (.4,.8)
    print(tri.area())
    print(tri.lengthSide(0,1))
    print(tri.lengthSide(2,1))
    print(tri.lengthSide(0,2))

    tri.drawVerts(win)
    tri.drawEdges(win)

    win.getMouse() # Pause to view result
    win.close()    # Close window when done


class Cube():

    def __init__(self, p1, p2, p3):
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3

    def volume(self):
        volume = abs(self.p2[0] - self.p1[0])
        volume *= abs(self.p3[1] - self.p2[1])
        volume *= abs(self.p3[2] - self.p1[2])
        print("Volume is " + str(volume))

def distance2D(p1, p2):
    return pow(pow(p2[0] - p1[0], 2) + pow(p2[1] - p1[1], 2), .5)

class Line():
    def __init__(self, slope, intercept):
        self.slope = slope
        self.intercept = intercept

    def getSlope(self):
        return self.slope

    def getSlope(p1, p2):
        # static function that returns the slope of the line between two points
        return -(p2[1] - p1[1]) / (p2[0] - p1[0])

    def intersection(self, other):
        if self.slope == other.slope:
            raise Exception("Lines parallel")
        x = (other.intercept - self.intercept)/(self.slope - other.slope)
        y = x*self.slope + self.intercept
        return (x,y)

class Tri():

    def __init__(self, verts):
        if len(verts) != 3:
            raise Exception("Invalid number of verts")
        self.verts = verts

    def lengthSide(self, vert1, vert2):
        return distance2D(self.verts[vert1], self.verts[vert2])

    def area(self):
        # calculate area using Heron's formula: http://jwilson.coe.uga.edu/emt725/Heron/HeronProofAlg.html
        a = self.lengthSide(0,1)
        b = self.lengthSide(1,2)
        c = self.lengthSide(2,0)
        s = (a + b + c) / 2
        area = pow(s*(s-a)*(s-b)*(s-c), 0.5)
        return area

    def drawVerts(self, win):
        p = Point(self.verts[0][0], self.verts[0][1])
        p.draw(win)
        p = Point(self.verts[1][0], self.verts[1][1])
        p.draw(win)
        p = Point(self.verts[2][0], self.verts[2][1])
        p.draw(win)

    def drawEdges(self, win):
        slope = Line.getSlope(self.verts[1], self.verts[2])
        for x in range(0, int(self.lengthSide(1,2))):
            p = Point(self.verts[1][0] + x*slope, self.verts[1][1] + x*slope)
            print(x)
            p.draw(win)


    #def inside(self, point):
main()
