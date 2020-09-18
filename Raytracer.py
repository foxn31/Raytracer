# Raytracer.py
# Written by Nick Fox using https://mitchellkember.com/blog/post/ray-tracer/
# as a reference

import os
from PIL import Image

class Vector():
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    # returns the mathematical magnitude of the vector
    def magnitude(self):
        return pow(self.x*self.x + self.y*self.y + self.z*self.z, 0.5)

    # returns a vector having the same direction as this one,
    # but with a magnitude of 1
    def normalize(self):
        magnitude = self.magnitude()
        return Vector(self.x/magnitude, self.y/magnitude, self.z/magnitude)

    # operator overloads and other useful functions:
    # add, sub, dot, cross, scalar multiply, direction, and distance
    def __add__(self, other_vector):
        return Vector(self.x + other_vector.x, self.y + other_vector.y,
            self.z + other_vector.z)

    def __sub__(self, other_vector):
        return Vector(self.x - other_vector.x, self.y - other_vector.y,
            self.z - other_vector.z)

    def dot(self, other_vector):
        return self.x*other_vector.x + self.y*other_vector.y + self.z*other_vector.z

    def cross(self, other_vector):
        x = self.y*other_vector.z - self.z*other_vector.y
        y = self.z*other_vector.x - self.x*other_vector.z
        z = self.x*other_vector.y - self.y*other_vector.x
        return Vector(x, y, z)

    def __mul__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return Vector(self.x*other, self.y*other, self.z*other)

    # returns the direction of a line passsing through two points in vector form
    def direction(v1, v2):
        direction = Vector(v2.x - v1.x, v2.y - v1.y, v2.z - v1.z)
        # normalize the direction before returning
        return direction.normalize()

    # returns the distance between two points in space, represented as vectors
    def distance(v1, v2):
        # distance formula, with minimal use of pow() for efficiency
        distance = (v2.x - v1.x)*(v2.x - v1.x)
        distance += (v2.y - v1.y)*(v2.y - v1.y)
        distance += (v2.z - v1.z)*(v2.z - v1.z)
        distance = pow(distance, 0.5)
        return distance

    # string representation tells x, y, and z components
    def __str__(self):
        return "(" + str(self.x) + "," + str(self.y) + "," + str(self.z) + ")"

class Ray():
    def __init__(self, pos, dir):
        self.pos = pos
        self.dir = dir

    # string representation tells position, direction, and magnitude
    def __str__(self):
        string = "p: " + str(self.pos) + "\n" + "d: " + str(self.dir)
        return string

class Sphere():
    def __init__(self, center, radius, material):
        self.center = center
        self.radius = radius
        self.material = material

    # determines the magnitude of the ray at which it intersects the Sphere
    # formula from mitchellkember, returns 0, 1, or 2 solutions
    def intersection(self, ray):
        # to solve, first calculate the terms of the quadratic equation to be solved
        a = 1
        b = (ray.dir*2).dot(ray.pos - self.center)
        c = (ray.pos - self.center).magnitude()*(ray.pos - self.center).magnitude()
        c -= self.radius*self.radius
        # perform quadratic formula on terms a, b, and c
        solutions = quadraticFormula(a, b, c)
        # since ray only points in one direction, remove any negative solutions
        validsolutions = []
        for solution in solutions:
            if solution >= 0:
                validsolutions.append(solution)
        # no intersection between ray and object
        if len(validsolutions) == 0:
            return None
        # return all points of intersection. the renderer will decide which one to use
        return validsolutions

    # given a point on the surface, return the normal vector
    # for a sphere, the point normal points away from the center
    def get_normal(self, point):
        return Vector.direction(point, self.center)

class Plane():
    def __init__(self, normal, scalar, material):
        self.normal = normal
        self.scalar = scalar
        self.material = material

    def intersection(self, ray):
        num = self.scalar - ray.pos.dot(self.normal)
        denom = ray.dir.dot(self.normal)
        # don't divide by zero (ray parallel to plane)
        # and don't return a negative value
        if denom == 0 or num / denom < 0:
            return None
        else:
            return [num / denom]

    # given a point on the surface, return the normal vector
    # for a plane, the normal is always the same
    def get_normal(self, point):
        return self.normal

# materials tell the renderer how to display an object
# albedo is an RGB vector telling the reflectivity for each color
class Material():
    def __init__(self, albedo):
        self.albedo = albedo

class Camera():
    def __init__(self, pos, dir, focal_length, width, height, pixel_width = 0.001):
        self.pos = pos
        self.dir = dir
        self.focal_length = focal_length
        self.width = width
        self.height = height
        self.pixel_width = pixel_width
    up = Vector(0, 0, 1)

class Light():
    def __init__(self, pos, power):
        self.pos = pos
        self.power = power

# solves a quadratic equation, returning 0, 1, or 2 solutions in a list
def quadraticFormula(a, b, c):
    rootTerm = b*b - 4*a*c
    # check for complex roots
    if rootTerm < 0:
        return []
    # solve for the two roots
    x1 = (-b + pow(rootTerm, 0.5)) / 2*a
    x2 = (-b - pow(rootTerm, 0.5)) / 2*a
    # check if the two roots are the same, in which case return just one
    if x1 == x2:
        return [x1]
    # return the roots as a tuple
    return [x1, x2]

# perform gamma correciton using sRGB color space (see mitchellkember)
def sRGB_gamma_correction(luminance):
    if (luminance <= 0.0031308):
        return 12.92*luminance
    else:
        return 1.055*pow(luminance, 1/2.4) - 0.055

# draws an image using ASCII text
def draw(image, width, height):
    for y in range(height):
        line = ""
        for x in range(width):
            if (image[x + width*y] == 0):
                line += "00"
            else:
                line += "11"
                #line += "▮"
        print(line)
    return

# draws an image using ASCII text and saves the result to a text file
def draw_to_file(image, width, height):
    with open("image.txt", "w+") as f:
        for y in range(height):
            line = ""
            for x in range(width):
                if image[x + width*y] == 0:
                    line += "77"
                elif image[x + width*y] == 1:
                    line += "  "
                    #line += "▮"
                elif image[x + width*y] == 2:
                    line += "88"
                elif image[x + width*y] == 3:
                    line += "//"
                elif image[x + width*y] == 4:
                    line += "**"
                elif image[x + width*y] == 5:
                    line += "--"
                elif image[x + width*y] == 10:
                    line += ".."
                else:
                    line += "55"
            f.write(line + "\n")
    return

def map_to_index(index, width):
    row = (index // width) // 2
    col = (index % width) // 2
    return int(row*width + col)

# perform post-processing techniques: exposure control, gamma correction,
# and anti-aliasing
def post_processing(sequence, width, height):
    # create a new list to hold the processed image
    exposure_gamma_sequence = []
    for pixel in sequence:
        # divide by 10 to control exposure
        r = pixel[0] / 10
        g = pixel[1] / 10
        b = pixel[2] / 10
        # gamma correction uses sRGB color space (see mitchellkember)
        r = sRGB_gamma_correction(r)
        g = sRGB_gamma_correction(g)
        b = sRGB_gamma_correction(b)
        # convert from [0.0, 1.0] to [0, 255] and round
        r = int(r * 255)
        g = int(g * 255)
        b = int(b * 255)
        exposure_gamma_sequence.append((r,g,b))

    return exposure_gamma_sequence

    old_index = 0
    new_width = width*2
    new_height = height*2

    interpolated_sequence = []
    # anti-aliasing using bilinear interpolation to resample
    for index in range(4 * len(exposure_gamma_sequence)):
        if (index // new_width) % 2 == 0:    # even-numbered row
            if index % 2 == 0:  # even-numbered column
                # the pixel is the same as the one in the original image
                #mapped_index = map_to_index(index, width)
                #print("index: " + str(index) + "  mapped_index: " + str(mapped_index))
                #interpolated_sequence.append(exposure_gamma_sequence[mapped_index])
                interpolated_sequence.append(exposure_gamma_sequence[old_index])
                old_index += 1
                #print(old_index)
            else:   # odd-numbered column
                # the pixel is the average of the pixels on either side,
                # unless at right edge, in which case its value is left neighbor's
                if (index % new_width) == (new_width - 1):
                    #mapped_index = map_to_index(index - 1, width)
                    mapped_index = old_index - 2
                    interpolated_sequence.append(exposure_gamma_sequence[mapped_index])
                else:
                    r = (exposure_gamma_sequence[(old_index - 1)][0]
                    + exposure_gamma_sequence[old_index][0]) // 2
                    g = (exposure_gamma_sequence[(old_index - 1)][1]
                    + exposure_gamma_sequence[old_index][1]) // 2
                    b = (exposure_gamma_sequence[(old_index - 1)][2]
                    + exposure_gamma_sequence[old_index][2]) // 2
                    interpolated_sequence.append((r,g,b))
        else:   # odd-numbered row
            if index % 2 == 0:  # even-numbered column
                # the pixel is the average of the pixels above and below,
                # unless at bottom edge, in which case its value is above neighbor's
                if (index // new_width) == (new_height - 1):
                    #mapped_index = map_to_index(index - new_width, width)
                    mapped_index = old_index - 1 - width
                    interpolated_sequence.append(exposure_gamma_sequence[mapped_index])
                else:
                    r = (exposure_gamma_sequence[(old_index - 1 - width)][0]
                    + exposure_gamma_sequence[old_index - 1 + width][0]) // 2
                    g = (exposure_gamma_sequence[(old_index - 1 - width)][1]
                    + exposure_gamma_sequence[old_index - 1 + width][1]) // 2
                    b = (exposure_gamma_sequence[(old_index - 1 - width)][2]
                    + exposure_gamma_sequence[old_index - 1 + width][2]) // 2
                    interpolated_sequence.append((r,g,b))
            else:   # odd-numbered column
                mapped_index = map_to_index(index, width)
                interpolated_sequence.append((0,0,0))

    return interpolated_sequence

def anti_aliasing(sequence, width, height):
    new_image = []
    new_width = width*2
    new_height = height*2

    for index, pixel in enumerate(sequence):
        # top-left pixel
        new_image.append(pixel)

        if index % width < width - 1:
            # top-middle
            r = (pixel[0] + sequence[index + 1][0]) // 2
            g = (pixel[1] + sequence[index + 1][1]) // 2
            b = (pixel[2] + sequence[index + 1][2]) // 2
            new_image.append((r,g,b))
        else:
            # if at the right edge, just duplicate the current pixel
            # top-middle
            r = (pixel[0])
            g = (pixel[1])
            b = (pixel[2])
            new_image.append((r,g,b))

        if index // width < height - 1:
            # middle-left
            r = (pixel[0] + sequence[index + width][0]) // 2
            g = (pixel[1] + sequence[index + width][1]) // 2
            b = (pixel[2] + sequence[index + width][2]) // 2
            new_image.append((r,g,b))

            # middle-middle
            # TODO: bilinear interpolation
            new_image.append((0,0,0))

        else:
            # if at the bottom edge, just duplicate the above pixel
            # middle-left
            r = sequence[index - width][0]
            g = sequence[index - width][1]
            b = sequence[index - width][2]
            new_image.append((r,g,b))

            # middle-middle
            # linear interpolation? just insert black pixels for now
            new_image.append((0,0,0))

    return new_image

def render_image(objects, light, camera):
    image = []
    # determine the positions of each pixel
    # this line gives the center of the image
    center = camera.pos + camera.dir*camera.focal_length
    # find the direction of the image's horizontal
    horizontal = camera.dir.cross(Camera.up).normalize()
    # each pixel needs to be offset from the center of the image
    for y in range(camera.height):
        for x in range(camera.width):
            h_offset = camera.pixel_width*(-0.5 + (camera.width/2 - x))
            v_offset = camera.pixel_width*(-0.5 + (camera.height/2 - y))
            # combine the horizontal and vertical offsets to find the position
            pixel = center + horizontal*h_offset + Vector(0, 0, v_offset)
            ray = Ray(camera.pos, Vector.direction(camera.pos, pixel))
            # only display the closest object to the camera
            closest_intersect = 10000
            closest_object = None
            for object in objects:
                intersects = object.intersection(ray)
                if intersects is not None:
                    # take the closest from the list of intersections returned
                    intersect = min(intersects)
                    if intersect < closest_intersect:
                        closest_intersect = intersect
                        closest_object = object
            if closest_object is None:
                image.append((0,0,0))
            else:
                # :cast another ray to see if this point is illuminated:
                # find the location of the intersect in 3D space
                intersect_point = ray.pos + ray.dir*closest_intersect
                # find the direction of the light and create the ray
                light_dir = Vector.direction(intersect_point, light.pos)
                ray = Ray(intersect_point, light_dir)
                # check all other objects in scene for occlusion
                distance_to_light = Vector.distance(intersect_point, light.pos)
                in_shadow = False
                for object in objects:
                    intersects = object.intersection(ray)
                    # if the ray intersects another object
                    if intersects is not None:
                        # take the farthest of the intersections to make sure
                        # it's not just intersecting with itself
                        # also compare the intersect with distance to light
                        # to ensure objects behind light aren't casting a shadow
                        intersect = max(intersects)
                        if intersect > .0001 and intersect <= distance_to_light:
                            # point is in shadow
                            in_shadow = True
                            image.append((0,0,0))
                            break
                if not in_shadow:
                    # point is not occluded;
                    # calculate color using Lambert (see mitchellkember):
                    # color = reflectivity*light_intensity*light_dir(dot)surface_normal
                    light_distance = Vector.distance(ray.pos, light.pos)
                    light_intensity = light.power / (light_distance * light_distance)
                    surface_normal = closest_object.get_normal(intersect_point)
                    lighting = light_intensity*(light_dir.dot(surface_normal))
                    #color = reflectivity*light_intensity*(light_dir.dot(surface_normal))
                    # perform calculation for each channel RGB, taking the
                    # absolute value since light_dir points in the direction
                    # opposite to the one we need
                    r = abs(closest_object.material.x*lighting)
                    g = abs(closest_object.material.y*lighting)
                    b = abs(closest_object.material.z*lighting)
                    image.append((r,g,b))
    return image

def main():
    ray = Ray(Vector(0,0,0), Vector(1,0,0))
    #print(ray)
    #draw(10,10)
    #print(quadraticFormula(1,0,-1))
    #print(Vector(1,0,1).dot(Vector(1, 1, 0)))
    material1 = Vector(1.0, 0.5, 0.5)
    material2 = Vector(1.0, 0.5, 1.0)
    material3 = Vector(0.8, 0.5, 0.0)
    material4 = Vector(0.5, 0.1, 1.0)
    material5 = Vector(0.8, 0.8, 0.8)
    material6 = Vector(0.1, 0.5, 0.5)
    sphere = Sphere(Vector(3, 0, 0), .2, material2)
    plane = Plane(Vector(1, 0, 0).normalize(), 2.7, material1)
    plane2 = Plane(Vector(0, 0, 1).normalize(), -0.6, material5)
    plane3 = Plane(Vector(1, 3, 3).normalize(), 0.2, material6)
    sphere2 = Sphere(Vector(3.3, .4, -.3), .1, material3)
    sphere3 = Sphere(Vector(2.9, -.3, .3), .05, material4)
    #sphere4 = Sphere(Vector(4.3, 1.3, 1.3), 0.1, material)
    #spheres = []
    #for x in range(36):
    #    sphere_pos = Vector(2.8, -0.7 + (x%6)*0.25, 0.6 - (x//6)*0.25)
    #    spheres.append(Sphere(sphere_pos, 0.1, material))
    #print(ray.dir*intersection(ray, sphere) + ray.pos)

    #objects = spheres + [plane]
    objects = [sphere, sphere2, sphere3, plane, plane2, plane3]
    #light = Light(Vector(4, 2, 2), 10)
    light = Light(Vector(4, 1, 1), 10)
    camera = Camera(Vector(11, 0, 0), Vector(-1, 0, 0).normalize(), 2, 400, 400)

    sequence = render_image(objects, light, camera)
    sequence = post_processing(sequence, camera.width, camera.height)
    sequence = anti_aliasing(sequence, camera.width, camera.height)

    #draw(sequence, camera.width, camera.height)
    #draw_to_file(sequence, camera.width, camera.height)
    img = Image.new("RGB", (camera.width*2, camera.height*2), (0,0,0))
    img.putdata(sequence)
    img.save("image.png", "PNG")

    os.system("image.png")

main()
